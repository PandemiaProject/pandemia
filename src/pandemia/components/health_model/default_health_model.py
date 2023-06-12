"""Default health model"""

import logging
import copy
import os
import csv
import numpy as np
from collections import defaultdict
from numpy import genfromtxt

from ctypes import c_void_p, c_double, c_int64, cdll

from ..health_model import HealthModel

log = logging.getLogger("default_health_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultHealthModel(HealthModel):
    """Default model of agent health.

    Instead of the compartmental description of health often applied in epidemic models, the default
    health model describes an agent's health with five attributes. These are sigma immunity
    (protection against infection), rho immunity (protection against severe outcomes),
    infectiousness, disease and strain. Low-level C functions, wrapped in Python functions,
    implement changes to agent health.

    Parameters:
    -----------
    config : Config
        A Pandemia Config object. A sub-config of the full config, containing the data used to
        configure this component.
    clock : Clock
        A Pandemia Clock object, discretizing the day.
    scale_factor : float
        The scale factor, coming from the full config.
    """

    def __init__(self, config, scale_factor, clock):
        """Initialize component"""
        super().__init__(config, scale_factor)

        self.scale_factor = scale_factor

        self.clock = clock

        self.update_health = self.lib.update_health
        self.transmission  = self.lib.transmission
        self.infect        = self.lib.infect

        self.update_health.restype = None
        self.transmission.restype  = None
        self.infect.restype        = None

        self.beta                             = np.array(config['beta'], dtype=np.float64)

        self.number_of_strains                = config['number_of_strains']
        self.num_initial_infections_by_region = config['num_initial_infections_by_region_by_strain']
        self.sir_rescaling                    = config['sir_rescaling']
        self.sir_rescaling_int                = int(self.sir_rescaling == True)

        self.age_mixing_matrices_directory    = config['age_mixing_matrices']
        self.age_group_interval               = config['age_group_interval']

        self.facemask_transmission_multiplier = config['facemask_transmission_multiplier']
        self.mutation_matrix                  = config['mutation_matrix']
        self.location_typ_multipliers         = config['location_typ_multipliers']
        self.immunity_period_days             = config['immunity_period_days']

        # Pre-existing immunity
        self.preexisting_sigma_multiplier     = config['preexisting_sigma_multiplier']
        self.preexisting_rho_multiplier       = config['preexisting_rho_multiplier']
        self.country_data_sigma_immunity_fp   = config['country_data_sigma_immunity_fp']
        self.country_data_rho_immunity_fp     = config['country_data_rho_immunity_fp']
        self.preexisting_sigma_immunity = {}
        self.preexisting_rho_immunity = {}
        if self.country_data_sigma_immunity_fp is not None:
            with open(self.country_data_sigma_immunity_fp, newline='') as csvfile:
                next(csvfile)
                region_data = csv.reader(csvfile, delimiter=',')
                for row in region_data:
                    iso2 = str(row[1])
                    self.preexisting_sigma_immunity[iso2] =\
                        [float(row[3 + r]) for r in range(self.number_of_strains)]
        if self.country_data_rho_immunity_fp is not None:
            with open(self.country_data_rho_immunity_fp, newline='') as csvfile:
                next(csvfile)
                region_data = csv.reader(csvfile, delimiter=',')
                for row in region_data:
                    iso2 = str(row[1])
                    self.preexisting_rho_immunity[iso2] =\
                        [float(row[3 + r]) for r in range(self.number_of_strains)]

        # Health presets
        if config['auto_generate_presets']:
            self.number_of_strains = 1
            sir_beta          = config['sir_beta']
            sir_gamma_inverse = config['sir_gamma_inverse']
            sir_disease_level = config['sir_disease_level']
            sir_max_dist_days = config['sir_max_dist_days']
            self.health_presets_config = defaultdict(dict)
            self.preset_weights_by_age = defaultdict(dict)
            self._generate_presets_and_weights(sir_beta, sir_gamma_inverse, sir_disease_level,
                                               sir_max_dist_days, self.clock)
        else:
            self.health_presets_config = config['health_presets']
            self.preset_weights_by_age = config['preset_weights_by_age']

        self.age_groups = list(self.preset_weights_by_age.keys())
        self.age_groups.sort()
        self.number_of_age_groups = len(self.age_groups)

        self.preset_ids_dict   = {}
        self.preset_names      = list(self.health_presets_config.keys())
        self.number_of_presets = len(self.preset_names)

        self.immunity_length = (self.clock.simulation_length_days // self.immunity_period_days) + 1
        self.immunity_period_ticks = self.immunity_period_days * self.clock.ticks_in_day

        self.number_of_rho_immunity_outcomes = self._get_number_of_rho_immunity_outcomes()
        self.max_preset_length_immunity      = self._get_max_preset_length_immunity()
        self.max_preset_length_health        = self._get_max_preset_length_health()

        self.preset_rho_immunity_failure_values =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains,
                      self.max_preset_length_immunity,
                      self.number_of_rho_immunity_outcomes), dtype=np.float64)

        self.preset_rho_immunity_failure_partitions =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains,
                      self.max_preset_length_immunity), dtype=np.int64)

        self.preset_rho_immunity_failure_lengths =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains), dtype=np.int64)

        self.preset_sigma_immunity_failure_values =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains,
                      self.max_preset_length_immunity), dtype=np.float64)

        self.preset_sigma_immunity_failure_partitions =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains,
                      self.max_preset_length_immunity), dtype=np.int64)

        self.preset_sigma_immunity_failure_lengths =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains), dtype=np.int64)

        self.preset_infectiousness_values =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=np.float64)

        self.preset_infectiousness_partitions =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=np.int64)

        self.preset_infectiousness_lengths =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains), dtype=np.int64)

        self.preset_disease_values =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=np.float64)

        self.preset_disease_partitions =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=np.int64)

        self.preset_disease_lengths =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains), dtype=np.int64)

        self.preset_strain_values =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=np.int64)

        self.preset_strain_partitions =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=np.int64)

        self.preset_strain_lengths =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains), dtype=np.int64)

        # Create presets, thereby filling the above arrays, and validate them
        self._get_presets(self.clock.ticks_in_day)
        self._validate_presets()

        # Create array of weighted presets by age group
        self.preset_weights = np.zeros((self.number_of_age_groups, self.number_of_presets))
        for a in range(self.number_of_age_groups):
            for preset_name in self.preset_names:
                id = self.preset_ids_dict[preset_name]
                age_group = self.age_groups[a]
                self.preset_weights[a][id] = self.preset_weights_by_age[age_group][preset_name]
        self.preset_weights /= np.sum(self.preset_weights, axis=1)[:, np.newaxis]

        # Flatten arrays
        self.mutation_matrix = np.asarray(self.mutation_matrix, dtype=np.float64).flatten()
        self.preset_infectiousness_lengths =\
            self.preset_infectiousness_lengths.flatten()
        self.preset_infectiousness_values =\
            self.preset_infectiousness_values.flatten()
        self.preset_infectiousness_partitions =\
            self.preset_infectiousness_partitions.flatten()
        self.preset_disease_lengths =\
            self.preset_disease_lengths.flatten()
        self.preset_disease_values =\
            self.preset_disease_values.flatten()
        self.preset_disease_partitions =\
            self.preset_disease_partitions.flatten()
        self.preset_strain_lengths =\
            self.preset_strain_lengths.flatten()
        self.preset_strain_values =\
            self.preset_strain_values.flatten()
        self.preset_strain_partitions =\
            self.preset_strain_partitions.flatten()
        self.preset_rho_immunity_failure_partitions =\
            self.preset_rho_immunity_failure_partitions.flatten()
        self.preset_rho_immunity_failure_values =\
            self.preset_rho_immunity_failure_values.flatten()
        self.preset_rho_immunity_failure_lengths =\
            self.preset_rho_immunity_failure_lengths.flatten()
        self.preset_sigma_immunity_failure_partitions =\
            self.preset_sigma_immunity_failure_partitions.flatten()
        self.preset_sigma_immunity_failure_values =\
            self.preset_sigma_immunity_failure_values.flatten()
        self.preset_sigma_immunity_failure_lengths =\
            self.preset_sigma_immunity_failure_lengths.flatten()

    def vectorize_component(self, vector_region):
        """Initializes vector region numpy arrays related to this component."""

        number_of_agents                = vector_region.number_of_agents
        number_of_locations             = vector_region.number_of_locations
        number_of_strains               = self.number_of_strains
        number_of_rho_immunity_outcomes = self.number_of_rho_immunity_outcomes
        immunity_length                 = self.immunity_length
        max_preset_length_health        = self.max_preset_length_health

        vector_region.health_age_group                 = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.location_transmission_multiplier = np.ones((number_of_locations), dtype=np.float64)
        vector_region.infection_event                  = np.full((number_of_agents), -1, dtype=np.int64)
        vector_region.presets                          = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.requesting_immunity_update       = np.zeros((number_of_agents), dtype=np.int64)

        vector_region.current_rho_immunity_failure   = np.ones((number_of_agents,
                                                                number_of_strains,
                                                                number_of_rho_immunity_outcomes),
                                                                dtype=np.float64)
        vector_region.current_sigma_immunity_failure = np.ones((number_of_agents,
                                                                number_of_strains),
                                                                dtype=np.float64)

        vector_region.current_infectiousness = np.zeros((number_of_agents), dtype=np.float64)
        vector_region.current_disease        = np.zeros((number_of_agents), dtype=np.float64)
        vector_region.current_strain         = np.full((number_of_agents), -1, dtype=np.int64)

        vector_region.rho_immunity_failure_values   = np.ones((number_of_agents,
                                                               number_of_strains,
                                                               immunity_length,
                                                               number_of_rho_immunity_outcomes),
                                                               dtype=np.float64)
        vector_region.sigma_immunity_failure_values = np.ones((number_of_agents,
                                                               number_of_strains,
                                                               immunity_length),
                                                               dtype=np.float64)

        vector_region.infectiousness_values     = np.zeros((number_of_agents,
                                                            max_preset_length_health), dtype=np.float64)
        vector_region.infectiousness_partitions = np.zeros((number_of_agents,
                                                            max_preset_length_health), dtype=np.int64)
        vector_region.infectiousness_lengths    = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.infectiousness_indexes    = np.ones((number_of_agents), dtype=np.int64)
        vector_region.disease_values            = np.zeros((number_of_agents,
                                                            max_preset_length_health), dtype=np.float64)
        vector_region.disease_partitions        = np.zeros((number_of_agents,
                                                            max_preset_length_health), dtype=np.int64)
        vector_region.disease_lengths           = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.disease_indexes           = np.ones((number_of_agents), dtype=np.int64)
        vector_region.strain_values             = np.full((number_of_agents,
                                                           max_preset_length_health), -1, dtype=np.int64)
        vector_region.strain_partitions         = np.zeros((number_of_agents,
                                                            max_preset_length_health), dtype=np.int64)
        vector_region.strain_lengths            = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.strain_indexes            = np.ones((number_of_agents), dtype=np.int64)

        self._get_age_mixing_matrix(vector_region)

    def initial_conditions_1(self, vector_region):
        """Establishes initial conditions for default health model."""

        # Complete assignment of default health function values
        num_r_vals = self.number_of_rho_immunity_outcomes
        vector_region.infectiousness_lengths[:] = 1
        vector_region.infectiousness_partitions[:, 0] = -1
        vector_region.disease_lengths[:] = 1
        vector_region.disease_partitions[:, 0] = -1
        vector_region.strain_lengths[:] = 1
        vector_region.strain_partitions[:, 0] = -1
        vector_region.rho_immunity_failure_values[:, :, :, num_r_vals - 1] = 0.0

        # Assign levels of pre-existing immunity
        if self.country_data_sigma_immunity_fp is not None:
            if vector_region.name in self.preexisting_sigma_immunity:
                levels = np.array(self.preexisting_sigma_immunity[vector_region.name])
                levels[levels < 1.0] *= self.preexisting_sigma_multiplier
                vector_region.sigma_immunity_failure_values[:, :, :] =\
                    levels[np.newaxis, :, np.newaxis]
        if (self.country_data_rho_immunity_fp is not None) and (num_r_vals > 1):
            if vector_region.name in self.preexisting_rho_immunity:
                levels = np.array(self.preexisting_rho_immunity[vector_region.name])
                levels[levels < 1.0] *= self.preexisting_rho_multiplier
                vector_region.rho_immunity_failure_values[:, :, :, 0] =\
                    levels[np.newaxis, :, np.newaxis]
        self.update(vector_region, 0)

    def initial_conditions_2(self, vector_region):
        """Establishes initial conditions for default health model."""

        # Determine location transmission multipliers
        if (self.location_typ_multipliers is not None):
            location_typ = vector_region.location_typ_strings
            multipliers = np.array([self.location_typ_multipliers.get(loc_typ, 1.0)
                                    for loc_typ in location_typ], dtype=np.float64)
            vector_region.location_transmission_multiplier = multipliers

        # Assign health presets for each agent by age
        age_group_indices =\
            np.searchsorted(np.asarray(self.age_groups), vector_region.age, side='right') - 1
        weights = self.preset_weights[age_group_indices]
        cum_weights = np.cumsum(weights, axis=1)
        random_probs = vector_region.prng.prng_np.random((vector_region.number_of_agents, 1))
        preset_ids = (random_probs < cum_weights).argmax(axis=1)
        vector_region.presets = preset_ids

        # Flatten arrays
        vector_region.infectiousness_partitions =\
            vector_region.infectiousness_partitions.flatten()
        vector_region.infectiousness_values =\
            vector_region.infectiousness_values.flatten()
        vector_region.disease_partitions =\
            vector_region.disease_partitions.flatten()
        vector_region.disease_values =\
            vector_region.disease_values.flatten()
        vector_region.strain_partitions =\
            vector_region.strain_partitions.flatten()
        vector_region.strain_values =\
            vector_region.strain_values.flatten()
        vector_region.current_rho_immunity_failure =\
            vector_region.current_rho_immunity_failure.flatten()
        vector_region.rho_immunity_failure_values =\
            vector_region.rho_immunity_failure_values.flatten()
        vector_region.current_sigma_immunity_failure =\
            vector_region.current_sigma_immunity_failure.flatten()
        vector_region.sigma_immunity_failure_values =\
            vector_region.sigma_immunity_failure_values.flatten()
        vector_region.subpopulation_mixing_matrix =\
            vector_region.subpopulation_mixing_matrix.flatten()

        # Initial infections
        num_initial_infections_rescaled = [0 for _ in range(self.number_of_strains)]
        if isinstance(self.num_initial_infections_by_region, dict):
            if vector_region.name in self.num_initial_infections_by_region:
                for s in range(self.number_of_strains):
                    num_initial_infections_rescaled[s] = int(self.scale_factor *\
                        self.num_initial_infections_by_region[vector_region.name][s])
        elif isinstance(self.num_initial_infections_by_region, int):
            for s in range(self.number_of_strains):
                num_initial_infections_rescaled[s] =\
                    int(self.scale_factor * self.num_initial_infections_by_region)
        elif isinstance(self.num_initial_infections_by_region, float):
            for s in range(self.number_of_strains):
                num_initial_infections_rescaled[s] =\
                    int(vector_region.number_of_agents * self.num_initial_infections_by_region)
        else:
            for s in range(self.number_of_strains):
                num_initial_infections_rescaled[s] = 1
        uninfected_population = np.where(vector_region.current_strain == -1)[0]
        for s in range(self.number_of_strains):
            num_initial_infections =\
                min(num_initial_infections_rescaled[s], vector_region.number_of_agents)
            initial_infections = vector_region.prng.prng_np.choice(uninfected_population,
                                                              num_initial_infections, replace=False)
            vector_region.infection_event[initial_infections] = s
            uninfected_population = np.setdiff1d(uninfected_population, initial_infections)

        self.infect_wrapper(vector_region, 0)

        self.update(vector_region, 0)

    def update(self, vector_region, t):
        """Updates agent health."""

        self.update_health(
            c_int64(t),
            c_int64(vector_region.number_of_agents),
            c_int64(self.number_of_strains),
            c_int64(self.immunity_length),
            c_int64(self.number_of_rho_immunity_outcomes),
            c_int64(self.max_preset_length_health),
            c_int64(self.immunity_period_ticks),
            c_void_p(vector_region.infectiousness_partitions.ctypes.data),
            c_void_p(vector_region.infectiousness_indexes.ctypes.data),
            c_void_p(vector_region.infectiousness_values.ctypes.data),
            c_void_p(vector_region.infectiousness_lengths.ctypes.data),
            c_void_p(vector_region.current_infectiousness.ctypes.data),
            c_void_p(vector_region.disease_partitions.ctypes.data),
            c_void_p(vector_region.disease_indexes.ctypes.data),
            c_void_p(vector_region.disease_values.ctypes.data),
            c_void_p(vector_region.disease_lengths.ctypes.data),
            c_void_p(vector_region.current_disease.ctypes.data),
            c_void_p(vector_region.strain_partitions.ctypes.data),
            c_void_p(vector_region.strain_indexes.ctypes.data),
            c_void_p(vector_region.strain_values.ctypes.data),
            c_void_p(vector_region.strain_lengths.ctypes.data),
            c_void_p(vector_region.current_strain.ctypes.data),
            c_void_p(vector_region.rho_immunity_failure_values.ctypes.data),
            c_void_p(vector_region.current_rho_immunity_failure.ctypes.data),
            c_void_p(vector_region.sigma_immunity_failure_values.ctypes.data),
            c_void_p(vector_region.current_sigma_immunity_failure.ctypes.data),
            c_void_p(vector_region.requesting_immunity_update.ctypes.data)
        )

    def dynamics(self, vector_region, t):
        """The transmission dynamics."""

        self.transmission(
            c_int64(self.number_of_strains),
            c_int64(vector_region.number_of_locations),
            c_int64(vector_region.number_of_agents),
            c_int64(vector_region.id),
            c_int64(self.sir_rescaling_int),
            c_int64(self.clock.ticks_in_day),
            c_int64(vector_region.number_of_subpopulations),
            c_void_p(self.beta.ctypes.data),
            c_void_p(vector_region.subpopulation_index.ctypes.data),
            c_void_p(vector_region.subpopulation_mixing_matrix.ctypes.data),
            c_double(self.facemask_transmission_multiplier),
            c_double(vector_region.current_region_transmission_multiplier),
            c_void_p(vector_region.current_strain.ctypes.data),
            c_void_p(vector_region.current_disease.ctypes.data),
            c_void_p(vector_region.current_facemask.ctypes.data),
            c_void_p(vector_region.current_location.ctypes.data),
            c_void_p(vector_region.current_region.ctypes.data),
            c_void_p(vector_region.current_infectiousness.ctypes.data),
            c_void_p(vector_region.location_transmission_multiplier.ctypes.data),
            c_void_p(self.mutation_matrix.ctypes.data),
            c_void_p(vector_region.current_sigma_immunity_failure.ctypes.data),
            c_void_p(vector_region.infection_event.ctypes.data),
            c_void_p(vector_region.random_state.ctypes.data)
        )

        self.infect_wrapper(vector_region, t)

    def infect_wrapper(self, vector_region, t):
        """The infection system."""

        self.infect(
            c_int64(t),
            c_int64(vector_region.number_of_agents),
            c_int64(self.number_of_rho_immunity_outcomes),
            c_int64(self.number_of_strains),
            c_int64(self.max_preset_length_health),
            c_int64(self.max_preset_length_immunity),
            c_int64(self.immunity_length),
            c_int64(self.immunity_period_ticks),
            c_void_p(vector_region.current_rho_immunity_failure.ctypes.data),
            c_void_p(vector_region.infection_event.ctypes.data),
            c_void_p(vector_region.infectiousness_partitions.ctypes.data),
            c_void_p(vector_region.infectiousness_values.ctypes.data),
            c_void_p(vector_region.infectiousness_lengths.ctypes.data),
            c_void_p(vector_region.infectiousness_indexes.ctypes.data),
            c_void_p(vector_region.disease_partitions.ctypes.data),
            c_void_p(vector_region.disease_values.ctypes.data),
            c_void_p(vector_region.disease_lengths.ctypes.data),
            c_void_p(vector_region.disease_indexes.ctypes.data),
            c_void_p(vector_region.strain_partitions.ctypes.data),
            c_void_p(vector_region.strain_values.ctypes.data),
            c_void_p(vector_region.strain_lengths.ctypes.data),
            c_void_p(vector_region.strain_indexes.ctypes.data),
            c_void_p(vector_region.rho_immunity_failure_values.ctypes.data),
            c_void_p(vector_region.sigma_immunity_failure_values.ctypes.data),
            c_void_p(vector_region.requesting_immunity_update.ctypes.data),
            c_void_p(self.preset_infectiousness_partitions.ctypes.data),
            c_void_p(self.preset_infectiousness_values.ctypes.data),
            c_void_p(self.preset_infectiousness_lengths.ctypes.data),
            c_void_p(self.preset_disease_partitions.ctypes.data),
            c_void_p(self.preset_disease_values.ctypes.data),
            c_void_p(self.preset_disease_lengths.ctypes.data),
            c_void_p(self.preset_strain_partitions.ctypes.data),
            c_void_p(self.preset_strain_values.ctypes.data),
            c_void_p(self.preset_strain_lengths.ctypes.data),
            c_void_p(self.preset_rho_immunity_failure_partitions.ctypes.data),
            c_void_p(self.preset_rho_immunity_failure_values.ctypes.data),
            c_void_p(self.preset_rho_immunity_failure_lengths.ctypes.data),
            c_void_p(self.preset_sigma_immunity_failure_partitions.ctypes.data),
            c_void_p(self.preset_sigma_immunity_failure_values.ctypes.data),
            c_void_p(self.preset_sigma_immunity_failure_lengths.ctypes.data),
            c_void_p(vector_region.presets.ctypes.data),
            c_void_p(vector_region.random_state.ctypes.data)
        )

    def _get_number_of_rho_immunity_outcomes(self):
        """Determines number of rho immunity outcomes in config"""

        num_rho_outcomes = len(self.health_presets_config[self.preset_names[0]])
        for p in self.health_presets_config:
            assert len(self.health_presets_config[p]) == num_rho_outcomes

        for preset_name in self.preset_names:
            for r in range(num_rho_outcomes):
                for s1 in range(self.number_of_strains):
                    for s2 in range(self.number_of_strains):
                        f = 'rho_immunity_failure'
                        items = len(self.health_presets_config[preset_name][r][f][s1][s2][1])
                        for i in range(items):
                            assert len(self.health_presets_config[preset_name][r][f][s1][s2][1][i])\
                                   == num_rho_outcomes

        return num_rho_outcomes

    def _get_max_preset_length_immunity(self):
        """Determine maximum length of immunity outcomes in config"""

        max_preset_length_immunity = 0
        for p in self.health_presets_config:
            for r in self.health_presets_config[p]:
                for s1 in range(self.number_of_strains):
                    for s2 in range(self.number_of_strains):
                        for f in ['rho_immunity_failure', 'sigma_immunity_failure']:
                            partition_data = self.health_presets_config[p][r][f][s1][s2][0]
                            values_data = self.health_presets_config[p][r][f][s1][s2][1]
                            assert len(partition_data) == len(values_data)
                            if len(partition_data) > max_preset_length_immunity:
                                max_preset_length_immunity = len(partition_data)

        return max_preset_length_immunity

    def _get_max_preset_length_health(self):
        """Determine maximum length of health outcomes in config"""

        max_preset_length_health = 0
        for p in self.health_presets_config:
            for r in self.health_presets_config[p]:
                for s in range(self.number_of_strains):
                    for f in ['infectiousness', 'disease', 'strain']:
                        partition_data = self.health_presets_config[p][r][f][s][0]
                        values_data = self.health_presets_config[p][r][f][s][1]
                        assert len(partition_data) == len(values_data)
                        if len(partition_data) > max_preset_length_health:
                            max_preset_length_health = len(partition_data)

        return max_preset_length_health

    def _get_presets(self, ticks_in_day):
        """Create health presets"""

        # Build infection response presets
        for p in self.health_presets_config:
            id = self.preset_names.index(p)
            self.preset_ids_dict[p] = id
            for r1 in self.health_presets_config[p]:
                f = 'rho_immunity_failure'
                for s1 in range(self.number_of_strains):
                    for s2 in range(self.number_of_strains):
                        partition_data = self.health_presets_config[p][r1][f][s1][s2][0]
                        values_data = self.health_presets_config[p][r1][f][s1][s2][1]
                        partition = [int(part * ticks_in_day) for part in partition_data]
                        values = [[float(val1) for val1 in val2] for val2 in values_data]
                        length = len(partition)
                        self.preset_rho_immunity_failure_lengths[id][r1][s1][s2] = length
                        for i in range(length):
                            self.preset_rho_immunity_failure_partitions[id][r1][s1][s2][i] =\
                                partition[i]
                            for r2 in range(self.number_of_rho_immunity_outcomes):
                                self.preset_rho_immunity_failure_values[id][r1][s1][s2][i][r2] =\
                                    values[i][r2]
                        for i in range(length, self.max_preset_length_immunity):
                            self.preset_rho_immunity_failure_partitions[id][r1][s1][s2][i] = -1
                            for r2 in range(self.number_of_rho_immunity_outcomes):
                                self.preset_rho_immunity_failure_values[id][r1][s1][s2][i][r2] = 1.0
                f = 'sigma_immunity_failure'
                for s1 in range(self.number_of_strains):
                    for s2 in range(self.number_of_strains):
                        partition_data = self.health_presets_config[p][r1][f][s1][s2][0]
                        values_data = self.health_presets_config[p][r1][f][s1][s2][1]
                        partition = [int(part * ticks_in_day) for part in partition_data]
                        values = [float(val) for val in values_data]
                        length = len(partition)
                        self.preset_sigma_immunity_failure_lengths[id][r1][s1][s2] = length
                        for i in range(length):
                            self.preset_sigma_immunity_failure_partitions[id][r1][s1][s2][i] =\
                                partition[i]
                            self.preset_sigma_immunity_failure_values[id][r1][s1][s2][i] =\
                                values[i]
                        for i in range(length, self.max_preset_length_immunity):
                            self.preset_sigma_immunity_failure_partitions[id][r1][s1][s2][i] = -1
                            self.preset_sigma_immunity_failure_values[id][r1][s1][s2][i] = 1.0
                f = 'infectiousness'
                for s in range(self.number_of_strains):
                    partition_data = self.health_presets_config[p][r1][f][s][0]
                    values_data = self.health_presets_config[p][r1][f][s][1]
                    partition = [int(part * ticks_in_day) for part in partition_data]
                    values = [float(val) for val in values_data]
                    length = len(partition)
                    self.preset_infectiousness_lengths[id][r1][s] = length
                    for i in range(length):
                        self.preset_infectiousness_partitions[id][r1][s][i] = partition[i]
                        self.preset_infectiousness_values[id][r1][s][i] = values[i]
                    for i in range(length, self.max_preset_length_health):
                        self.preset_infectiousness_partitions[id][r1][s][i] = -1
                        self.preset_infectiousness_values[id][r1][s][i] = 0.0
                f = 'disease'
                for s in range(self.number_of_strains):
                    partition_data = self.health_presets_config[p][r1][f][s][0]
                    values_data = self.health_presets_config[p][r1][f][s][1]
                    partition = [int(part * ticks_in_day) for part in partition_data]
                    values = [float(val) for val in values_data]
                    length = len(partition)
                    self.preset_disease_lengths[id][r1][s] = length
                    for i in range(length):
                        self.preset_disease_partitions[id][r1][s][i] = partition[i]
                        self.preset_disease_values[id][r1][s][i] = values[i]
                    for i in range(length, self.max_preset_length_health):
                        self.preset_disease_partitions[id][r1][s][i] = -1
                        self.preset_disease_values[id][r1][s][i] = 0.0
                f = 'strain'
                for s in range(self.number_of_strains):
                    partition_data = self.health_presets_config[p][r1][f][s][0]
                    values_data = self.health_presets_config[p][r1][f][s][1]
                    partition = [int(part * ticks_in_day) for part in partition_data]
                    values = [float(val) for val in values_data]
                    length = len(partition)
                    self.preset_strain_lengths[id][r1][s] = length
                    for i in range(length):
                        self.preset_strain_partitions[id][r1][s][i] = partition[i]
                        self.preset_strain_values[id][r1][s][i] = values[i]
                    for i in range(length, self.max_preset_length_health):
                        self.preset_strain_partitions[id][r1][s][i] = -1
                        self.preset_strain_values[id][r1][s][i] = -1

    def _validate_presets(self):
        """Validate presets"""

        for id in range(self.number_of_presets):
            for r in range(self.number_of_rho_immunity_outcomes):
                for s in range(self.number_of_strains):
                    infectiousness_length = self.preset_infectiousness_lengths[id][r][s]
                    disease_length = self.preset_disease_lengths[id][r][s]
                    strain_length = self.preset_strain_lengths[id][r][s]
                    final_infectiousness_partition =\
                        self.preset_infectiousness_partitions[id][r][s][infectiousness_length - 1]
                    final_disease_partition =\
                        self.preset_disease_partitions[id][r][s][disease_length - 1]
                    final_strain_partition =\
                        self.preset_strain_partitions[id][r][s][strain_length - 1]
                    assert final_infectiousness_partition <= final_strain_partition
                    assert final_disease_partition <= final_strain_partition

    def _get_age_mixing_matrix(self, vector_region):
        """Extracts age mixing matrix from data (optional)"""

        if (self.age_mixing_matrices_directory is not None) and self.sir_rescaling:
            for root, dir, files in os.walk(self.age_mixing_matrices_directory):
                if vector_region.name + ".csv" in files:
                    mat = genfromtxt(self.age_mixing_matrices_directory +\
                                     vector_region.name + ".csv", skip_header = 1,
                                     delimiter=',', dtype=np.float64)
                    assert mat.shape[0] == mat.shape[1]
                    vector_region.number_of_subpopulations = mat.shape[0]

                    # Symmetrize matrix
                    mat = (mat + np.transpose(mat)) / 2

                    # Normalize rows to probability vectors
                    row_sums = mat.sum(axis=1)
                    mat = mat / row_sums[:, np.newaxis]

                    # Scale rows to take account of the fact that some age groups have more
                    # contacts in total than other age groups
                    row_sums = row_sums / np.amax(row_sums)
                    mat = mat * row_sums[:, np.newaxis]

                    # Calculate share of population in each age mixing group
                    share_in_age_group = np.zeros((mat.shape[0]), dtype=np.float64)
                    for n in range(vector_region.number_of_agents):
                        age_group = min(vector_region.age[n] // self.age_group_interval,
                                        mat.shape[0] - 1)
                        share_in_age_group[age_group] += 1
                    share_in_age_group /= share_in_age_group.sum(axis=0)

                    # Adjust final matrix so that the inner product of the row sums and the vector
                    # share_in_age_group is equal to 1.0, that is, the same for each region
                    dim = vector_region.number_of_subpopulations
                    epsilon = (1.0 - np.inner(mat.sum(axis=1), share_in_age_group)) / dim
                    mat += np.full((dim, dim), epsilon, dtype=np.float64)
                    vector_region.subpopulation_mixing_matrix = mat

                else:
                    vector_region.subpopulation_mixing_matrix = np.ones((1, 1), dtype=np.float64)
                    vector_region.number_of_subpopulations = 1
            vector_region.age_group = np.zeros((vector_region.number_of_agents), dtype=np.int64)
            for n in range(vector_region.number_of_agents):
                age_group = min(vector_region.age[n] // self.age_group_interval,
                                vector_region.number_of_subpopulations - 1)
                vector_region.subpopulation_index[n] = age_group
        else:
            vector_region.subpopulation_mixing_matrix = np.ones((1, 1), dtype=np.float64)
            vector_region.subpopulation_index =\
                np.zeros((vector_region.number_of_agents), dtype=np.int64)
            vector_region.number_of_subpopulations = 1

    def _generate_presets_and_weights(self, sir_beta, sir_gamma_inverse,
                                      sir_disease_level, sir_max_dist_days, clock):
        """Generates simple health presets using geometric distribution"""

        ticks_in_day = clock.ticks_in_day
        max_day = clock.simulation_length_days

        for t in range(1, sir_max_dist_days * ticks_in_day):
            preset_name = "preset_" + str(t - 1)
            self.health_presets_config[preset_name] =\
                {
                    0:
                    {
                        "rho_immunity_failure":
                            [[[[-1, t / ticks_in_day, max_day], [[0.0], [0.0], [0.0]]]]],
                        "sigma_immunity_failure":
                            [[[[-1, t / ticks_in_day, max_day], [1.0, 0.0, 1.0]]]],
                        "infectiousness":
                            [[[-1, 0, t / ticks_in_day], [0.0, sir_beta, 0.0]]],
                        "disease":
                            [[[-1, 0, t / ticks_in_day], [0.0, sir_disease_level, 0.0]]],
                        "strain":
                            [[[-1, 0, t / ticks_in_day], [-1, 0, -1]]]
                    }
                }
            self.preset_weights_by_age[0][preset_name] =\
                ((1 - (1 / (sir_gamma_inverse * ticks_in_day)))**(t - 1)) *\
                    (1 / (sir_gamma_inverse * ticks_in_day))
