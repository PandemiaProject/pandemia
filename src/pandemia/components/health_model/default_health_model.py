"""Default health model"""

import logging
import copy
import os
import csv
import numpy as np
from collections import defaultdict
from numpy import genfromtxt

from ctypes import c_void_p, c_double, c_int, cdll

from pandemia.components.health_model import HealthModel

log = logging.getLogger("default_health_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultHealthModel(HealthModel):
    """Default model of agent health. This model abandons the compartmental description of health
    often applied in epidemic models. Instead, an agent's health is determined by five parameters:
    
        - sigma immunity
        - rho immunity
        - infectiousness
        - disease
        - strain
    
    """

    def __init__(self, config, scale_factor, clock):
        """Initialize component"""
        super().__init__(config, scale_factor)

        self.scale_factor = scale_factor

        self.clock = clock

        lib = cdll.LoadLibrary("./src/pandemia/components/health_model/"
                                "default_health_model_functions.dll")

        self.update_health = lib.update_health
        self.transmission  = lib.transmission
        self.infect        = lib.infect

        self.update_health.restype = None
        self.transmission.restype  = None
        self.infect.restype        = None

        self.beta                             = np.array(config['beta'], dtype=float)

        self.number_of_strains                = config['number_of_strains']
        self.num_initial_infections_by_region = config['num_initial_infections_by_region_by_strain']

        self.sir_rescaling                    = config['sir_rescaling']
        self.sir_rescaling_int                = int(self.sir_rescaling == True)

        self.age_mixing_matrices_directory    = config['age_mixing_matrices']
        self.age_group_interval               = config['age_group_interval']

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

        self.facemask_transmission_multiplier = config['facemask_transmission_multiplier']
        self.mutation_matrix                  = config['mutation_matrix']
        self.location_typ_multipliers         = config['location_typ_multipliers']
        self.immunity_period_days             = config['immunity_period_days']

        self.age_groups = list(self.preset_weights_by_age.keys())
        self.age_groups.sort()

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
                      self.number_of_rho_immunity_outcomes), dtype=float)

        self.preset_rho_immunity_failure_partitions =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains,
                      self.max_preset_length_immunity), dtype=int)

        self.preset_rho_immunity_failure_lengths =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains), dtype=int)

        self.preset_sigma_immunity_failure_values =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains,
                      self.max_preset_length_immunity), dtype=float)

        self.preset_sigma_immunity_failure_partitions =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains,
                      self.max_preset_length_immunity), dtype=int)

        self.preset_sigma_immunity_failure_lengths =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.number_of_strains), dtype=int)

        self.preset_infectiousness_values =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=float)

        self.preset_infectiousness_partitions =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=int)

        self.preset_infectiousness_lengths =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains), dtype=int)

        self.preset_disease_values =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=float)

        self.preset_disease_partitions =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=int)

        self.preset_disease_lengths =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains), dtype=int)

        self.preset_strain_values =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=int)

        self.preset_strain_partitions =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains,
                      self.max_preset_length_health), dtype=int)

        self.preset_strain_lengths =\
            np.zeros((self.number_of_presets,
                      self.number_of_rho_immunity_outcomes,
                      self.number_of_strains), dtype=int)

        # Create presets, thereby filling the above arrays, and validate them
        self._get_presets(self.clock.ticks_in_day)
        self._validate_presets()

        # Flatten arrays
        self.mutation_matrix = np.asarray(self.mutation_matrix, dtype=float).flatten()
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
        """Initializes vector region numpy arrays related to this component"""

        number_of_agents                = vector_region.number_of_agents
        number_of_locations             = vector_region.number_of_locations
        number_of_strains               = self.number_of_strains
        number_of_rho_immunity_outcomes = self.number_of_rho_immunity_outcomes
        immunity_length                 = self.immunity_length
        max_preset_length_health        = self.max_preset_length_health

        vector_region.health_age_group                 = np.zeros((number_of_agents), dtype=int)
        vector_region.location_transmission_multiplier = np.ones((number_of_locations), dtype=float)
        vector_region.infection_event                  = np.full((number_of_agents), -1, dtype=int)
        vector_region.presets                          = np.zeros((number_of_agents), dtype=int)
        vector_region.requesting_immunity_update       = np.zeros((number_of_agents), dtype=int)

        vector_region.current_rho_immunity_failure   = np.ones((number_of_agents,
                                                                number_of_strains,
                                                                number_of_rho_immunity_outcomes),
                                                                dtype=float)
        vector_region.current_sigma_immunity_failure = np.ones((number_of_agents,
                                                                number_of_strains),
                                                                dtype=float)

        vector_region.current_infectiousness = np.zeros((number_of_agents), dtype=float)
        vector_region.current_disease        = np.zeros((number_of_agents), dtype=float)
        vector_region.current_strain         = np.full((number_of_agents), -1, dtype=int)

        vector_region.rho_immunity_failure_values   = np.ones((number_of_agents,
                                                               number_of_strains,
                                                               immunity_length,
                                                               number_of_rho_immunity_outcomes),
                                                               dtype=float)
        vector_region.sigma_immunity_failure_values = np.ones((number_of_agents,
                                                               number_of_strains,
                                                               immunity_length),
                                                               dtype=float)

        vector_region.infectiousness_values     = np.zeros((number_of_agents,
                                                            max_preset_length_health), dtype=float)
        vector_region.infectiousness_partitions = np.zeros((number_of_agents,
                                                            max_preset_length_health), dtype=int)
        vector_region.infectiousness_lengths    = np.zeros((number_of_agents), dtype=int)
        vector_region.infectiousness_indexes    = np.ones((number_of_agents), dtype=int)
        vector_region.disease_values            = np.zeros((number_of_agents,
                                                            max_preset_length_health), dtype=float)
        vector_region.disease_partitions        = np.zeros((number_of_agents,
                                                            max_preset_length_health), dtype=int)
        vector_region.disease_lengths           = np.zeros((number_of_agents), dtype=int)
        vector_region.disease_indexes           = np.ones((number_of_agents), dtype=int)
        vector_region.strain_values             = np.full((number_of_agents,
                                                           max_preset_length_health), -1, dtype=int)
        vector_region.strain_partitions         = np.zeros((number_of_agents,
                                                            max_preset_length_health), dtype=int)
        vector_region.strain_lengths            = np.zeros((number_of_agents), dtype=int)
        vector_region.strain_indexes            = np.ones((number_of_agents), dtype=int)

        self._get_age_mixing_matrix(vector_region)

    def initial_conditions(self, vector_region):
        """Updates health functions"""

        # Complete assignment of default health function values
        num_r_vals = self.number_of_rho_immunity_outcomes
        for n in range(vector_region.number_of_agents):
            vector_region.infectiousness_lengths[n] = 1
            vector_region.infectiousness_partitions[n][0] = -1
            vector_region.disease_lengths[n] = 1
            vector_region.disease_partitions[n][0] = -1
            vector_region.strain_lengths[n] = 1
            vector_region.strain_partitions[n][0] = -1
            for s in range(self.number_of_strains):
                for index in range(self.immunity_length):
                    vector_region.rho_immunity_failure_values[n][s][index][num_r_vals - 1] = 0.0

        # Assign levels of pre-existing immunity
        if self.country_data_sigma_immunity_fp is not None:
            if vector_region.name in self.preexisting_sigma_immunity:
                for s in range(self.number_of_strains):
                    level = self.preexisting_sigma_immunity[vector_region.name][s]
                    if level < 1.0:
                        level *= self.preexisting_sigma_multiplier
                        for n in range(vector_region.number_of_agents):
                            for index in range(self.immunity_length):
                                vector_region.sigma_immunity_failure_values[n][s][index] = level
        if (self.country_data_rho_immunity_fp is not None) and (num_r_vals > 1):
            if vector_region.name in self.preexisting_rho_immunity:
                for s in range(self.number_of_strains):
                    level = self.preexisting_rho_immunity[vector_region.name][s]
                    if level < 1.0:
                        level *= self.preexisting_rho_multiplier
                        for n in range(vector_region.number_of_agents):
                            for index in range(self.immunity_length):
                                vector_region.rho_immunity_failure_values[n][s][index][0] = level

        self.update(vector_region, 0)

        # Determine location transmission multipliers
        for m in range(vector_region.number_of_locations):
            location_typ = vector_region.location_typ_strings[m]
            if (self.location_typ_multipliers is not None) and\
                (location_typ in self.location_typ_multipliers):
                multiplier = self.location_typ_multipliers[location_typ]
            else:
                multiplier = 1.0
            vector_region.location_transmission_multiplier[m] = multiplier

        # Assign health presets for each agent by age
        for n in range(vector_region.number_of_agents):
            age = vector_region.age[n]
            age_group = self._determine_age_group(age, self.age_groups)
            preset_name =\
                self._assign_preset_name(age_group, vector_region.prng, vector_region.name)
            preset_id = self.preset_ids_dict[preset_name]
            vector_region.presets[n] = preset_id

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
        uninfected_population =\
            copy.deepcopy([n for n in range(vector_region.number_of_agents) if
                        vector_region.current_strain[n] == -1])
        for s in range(self.number_of_strains):
            num_initial_infections =\
                min(num_initial_infections_rescaled[s], vector_region.number_of_agents)
            initial_infections =\
                vector_region.prng.random_sample(uninfected_population,
                                                 num_initial_infections)
            for n in initial_infections:
                vector_region.infection_event[n] = s
                uninfected_population.remove(n)
     
        self.infect_wrapper(vector_region, 0)

        self.update(vector_region, 0)

    def update(self, vector_region, t):
        """Updates health functions"""

        self.update_health(
            c_int(t),
            c_int(vector_region.number_of_agents),
            c_int(self.number_of_strains),
            c_int(self.immunity_length),
            c_int(self.number_of_rho_immunity_outcomes),
            c_int(self.max_preset_length_health),
            c_int(self.immunity_period_ticks),
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
        """Changes to agent health"""

        self.transmission(
            c_int(self.number_of_strains),
            c_int(vector_region.number_of_locations),
            c_int(vector_region.number_of_agents),
            c_int(vector_region.id),
            c_int(self.sir_rescaling_int),
            c_int(self.clock.ticks_in_day),
            c_int(vector_region.number_of_subpopulations),
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
        """Changes to agent health"""

        self.infect(
            c_int(t),
            c_int(vector_region.number_of_agents),
            c_int(self.number_of_rho_immunity_outcomes),
            c_int(self.number_of_strains),
            c_int(self.max_preset_length_health),
            c_int(self.max_preset_length_immunity),
            c_int(self.immunity_length),
            c_int(self.immunity_period_ticks),
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

    def _assign_preset_name(self, age_group, prng, region_name):
        """Assigns agents infection response preset names by age"""

        preset_name = prng.multinoulli_dict(self.preset_weights_by_age[age_group])

        return preset_name

    def _determine_age_group(self, age, age_groups):
        """Determine to which age group the agent belongs"""

        if age >= age_groups[-1]:
            return age_groups[-1]
        elif age < age_groups[0]:
            return age_groups[0]
        else:
            i = 0
            while age < age_groups[i] or age >= age_groups[i+1]:
                i = i + 1
            return age_groups[i]

    def _get_age_mixing_matrix(self, vector_region):
        """Extracts age mixing matrix from data (optional)"""

        if (self.age_mixing_matrices_directory is not None) and self.sir_rescaling:
            for root, dir, files in os.walk(self.age_mixing_matrices_directory):
                if vector_region.name + ".csv" in files:
                    mat = genfromtxt(self.age_mixing_matrices_directory +\
                                     vector_region.name + ".csv", skip_header = 1,
                                     delimiter=',', dtype=float)
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
                    share_in_age_group = np.zeros((mat.shape[0]), dtype=float)
                    for n in range(vector_region.number_of_agents):
                        age_group = min(vector_region.age[n] // self.age_group_interval,
                                        mat.shape[0] - 1)
                        share_in_age_group[age_group] += 1
                    share_in_age_group /= share_in_age_group.sum(axis=0)

                    # Adjust final matrix so that the inner product of the row sums and the vector
                    # share_in_age_group is equal to 1.0, that is, the same for each region
                    dim = vector_region.number_of_subpopulations
                    epsilon = (1.0 - np.inner(mat.sum(axis=1), share_in_age_group)) / dim
                    mat += np.full((dim, dim), epsilon, dtype=float)
                    vector_region.subpopulation_mixing_matrix = mat

                else:
                    vector_region.subpopulation_mixing_matrix = np.ones((1, 1), dtype=float)
                    vector_region.number_of_subpopulations = 1
            vector_region.age_group = np.zeros((vector_region.number_of_agents), dtype=int)
            for n in range(vector_region.number_of_agents):
                age_group = min(vector_region.age[n] // self.age_group_interval,
                                vector_region.number_of_subpopulations - 1)
                vector_region.subpopulation_index[n] = age_group
        else:
            vector_region.subpopulation_mixing_matrix = np.ones((1, 1), dtype=float)
            vector_region.subpopulation_index =\
                np.zeros((vector_region.number_of_agents), dtype=int)
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
