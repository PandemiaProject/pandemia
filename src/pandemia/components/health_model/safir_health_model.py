"""Safir health model"""

import logging
import copy
import os
import csv
import numpy as np
from collections import defaultdict
from numpy import genfromtxt

from ctypes import c_void_p, c_double, c_int64, cdll

from ..health_model import HealthModel

log = logging.getLogger("safir_health_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class SafirHealthModel(HealthModel):
    """A version of the Safir health model.

    To see the Safir3 model, follow the link: https://github.com/mrc-ide/safir3.

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

        self.number_of_strains = 1
        self.immunity_length = 1
        self.immunity_period_ticks = 1
        self.mutation_matrix = [[1.0]]
        self.number_of_rho_immunity_outcomes = 1

        self.scale_factor = scale_factor

        self.clock = clock
        self.number_of_days = self.clock.simulation_length_days
        self.ticks_in_day = self.clock.ticks_in_day

        self.transmission = self.lib.safir3_transmission_model
        self.transmission.restype = None

        self.beta                             = np.array(config['beta'], dtype=np.float64)
        self.num_initial_infections_by_region = config['num_initial_infections_by_region']
        self.facemask_transmission_multiplier = config['facemask_transmission_multiplier']
        self.location_typ_multipliers         = config['location_typ_multipliers']

        self.sir_rescaling                    = config['sir_rescaling']
        self.sir_rescaling_int                = int(self.sir_rescaling == True)

        self.number_of_age_groups = 17

        self.age_mixing_matrices_directory    = config['age_mixing_matrices']
        self.age_group_interval               = config['age_group_interval']

        # Disease
        self.infectiousness_profile =\
            np.array([0.000395985, 0.001768802, 0.125221560, 0.222678871,
                      0.210222459, 0.176880156, 0.125221560, 0.070417258,
                      0.039598534, 0.017688016, 0.007041726, 0.002226789,
                      0.000994670], dtype=np.float64)
        self.infectious_period_days = self.infectiousness_profile.size
        self.infectious_period_ticks = self.infectious_period_days * self.ticks_in_day
        self.vfr = np.ones(self.number_of_days, dtype=np.float64)
        self.dr_vec = self._build_dr_vec()

    def vectorize_component(self, vector_region):
        """Initializes vector region numpy arrays related to this component"""

        number_of_agents = vector_region.number_of_agents

        vector_region.ab_titre =\
            np.full(number_of_agents, -np.Inf, dtype=np.float64)
        vector_region.infection_number =\
            np.zeros(number_of_agents, dtype=np.int64)
        vector_region.days_since_last_infection =\
            np.full(number_of_agents, self.number_of_days + 1, dtype=np.int64)
        vector_region.ef_infection =\
            np.ones(number_of_agents, dtype=np.float64)
        vector_region.ef_transmission =\
            np.ones(number_of_agents, dtype=np.float64)
        vector_region.ef_severe =\
            np.ones(number_of_agents, dtype=np.float64)
        vector_region.current_strain =\
            np.full(number_of_agents, -1, dtype=np.int64)
        vector_region.current_disease =\
            np.full(number_of_agents, 0, dtype=np.float64)
        vector_region.current_infectiousness =\
            np.full(number_of_agents, 0, dtype=np.float64)
        vector_region.start_t =\
            np.full(number_of_agents, -1, dtype=np.int64)
        vector_region.t_for_recovery =\
            np.full(number_of_agents, -1, dtype=np.int64)
        vector_region.infection_event =\
            np.full(number_of_agents, -1, dtype=np.int64)
        vector_region.location_transmission_multiplier =\
            np.ones(vector_region.number_of_locations, dtype=np.float64)

        self._load_contact_matrix(vector_region)

    def initial_conditions(self, vector_region):
        """Updates health functions"""

        number_of_agents = vector_region.number_of_agents

        # Determine location transmission multipliers
        for m in range(vector_region.number_of_locations):
            location_typ = vector_region.location_typ_strings[m]
            if (self.location_typ_multipliers is not None) and\
                (location_typ in self.location_typ_multipliers):
                multiplier = self.location_typ_multipliers[location_typ]
            else:
                multiplier = 1.0
            vector_region.location_transmission_multiplier[m] = multiplier

        # Number of initial infections
        num_initial_infections_rescaled = 0
        if isinstance(self.num_initial_infections_by_region, dict):
            if vector_region.name in self.num_initial_infections_by_region:
                num_initial_infections_rescaled = int(self.scale_factor *\
                    self.num_initial_infections_by_region[vector_region.name])
        elif isinstance(self.num_initial_infections_by_region, int):
            num_initial_infections_rescaled =\
                int(self.scale_factor * self.num_initial_infections_by_region)
        elif isinstance(self.num_initial_infections_by_region, float):
            num_initial_infections_rescaled =\
                int(number_of_agents * self.num_initial_infections_by_region)
        else:
            num_initial_infections_rescaled = 1
        num_initial_infections =\
            min(num_initial_infections_rescaled, number_of_agents)

        # Initial infections
        initial_infections =\
            vector_region.prng.random_sample(range(number_of_agents), num_initial_infections)
        vector_region.infection_event[initial_infections] = 0
        self.update(vector_region, 0)

    def update(self, vector_region, t):
        """Updates health functions"""

        tick = t % self.ticks_in_day
        day = t // self.ticks_in_day

        # Infections
        indexes = np.argwhere(vector_region.infection_event == 0).flatten()
        vector_region.t_for_recovery[indexes] = t + self.infectious_period_ticks
        vector_region.start_t[indexes] = t
        vector_region.current_strain[indexes] = 0
        vector_region.infection_event = np.full(vector_region.number_of_agents, -1, dtype=np.int64)

        # Recoveries
        indexes = np.argwhere(vector_region.t_for_recovery == t).flatten()
        self._recovery(vector_region, indexes, day)

        # Infectiousness
        vector_region.current_infectiousness =\
            np.full(vector_region.number_of_agents, 0, dtype=np.float64)
        indexes = np.argwhere(vector_region.current_strain == 0).flatten()
        vector_region.current_infectiousness[indexes] =\
            self.infectiousness_profile[(t - vector_region.start_t[indexes]) // self.ticks_in_day]

        # Daily updates
        if tick == 0:

            # Update immune efficacies
            vector_region.ef_infection =\
                self._efficacy_infection(vector_region.ab_titre, self.vfr, day)
            vector_region.ef_transmission =\
                self._efficacy_transmission(vector_region.ab_titre, self.vfr, day)
            vector_region.ef_severe =\
                self._efficacy_severe(vector_region.ab_titre, self.vfr, day,
                                      vector_region.ef_infection)

            # Update antibody titres
            self._update_ab_titre(vector_region, self.dr_vec, self.number_of_days)

            # Update days since last infection
            indexes = np.argwhere(vector_region.days_since_last_infection <
                                  self.number_of_days + 1).flatten()
            vector_region.days_since_last_infection[indexes] += 1

    def dynamics(self, vector_region, t):
        """Changes to agent health"""

        self.transmission(
            c_int64(vector_region.number_of_locations),
            c_int64(vector_region.number_of_agents),
            c_int64(vector_region.id),
            c_int64(self.sir_rescaling_int),
            c_int64(self.ticks_in_day),
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
            c_void_p(vector_region.ef_transmission.ctypes.data),
            c_void_p(vector_region.ef_infection.ctypes.data),
            c_void_p(vector_region.infection_event.ctypes.data),
            c_void_p(vector_region.random_state.ctypes.data)
        )

    def _build_dr_vec(self):
        """Builds dr_vec array"""

        period_s = 250
        t_period_l = 365
        hl_s = 108
        hl_l = 3650

        time_to_decay = t_period_l - period_s
        dr_s = -np.log(2)/hl_s
        dr_l = -np.log(2)/hl_l

        dr_vec = np.concatenate((np.full(period_s, dr_s, dtype=np.float64),
                                np.linspace(start=dr_s, stop=dr_l, num=time_to_decay),
                                np.array([dr_l], dtype=np.float64)), dtype=np.float64)

        return dr_vec

    def _efficacy_infection(self, ab_titre, vfr, day):
        """Converts ab_titre into efficacy against infection"""

        ab_50 = 0.2
        k = 2.94
        eps = np.finfo(float).eps

        ab_titre = np.exp(ab_titre) / vfr[day]
        ab_titre[ab_titre < eps] = eps
        ef_infection = 1 - (1 / (1 + np.exp(-k * (np.log10(ab_titre) - np.log10(ab_50)))))

        return ef_infection

    def _efficacy_transmission(self, ab_titre, vfr, day):
        """Converts ab_titre into efficacy against transmission"""

        ab_50 = 0.2
        k = 2.94
        nt_transmission_factor = 12
        eps = np.finfo(float).eps

        ab_titre = np.exp(ab_titre) / vfr[day]
        ab_titre[ab_titre < eps] = eps
        ef_transmission =\
            1 - (1 / (1 + np.exp(-k * (np.log10(ab_titre / nt_transmission_factor) - np.log10(ab_50)))))

        return ef_transmission

    def _efficacy_severe(self, ab_titre, vfr, day, ef_infection):
        """Converts ab_titre into efficacy against severe outcomes"""

        ab_50_severe = 0.03
        k = 2.94
        eps = np.finfo(float).eps

        ab_titre = np.exp(ab_titre) / vfr[day]
        ab_titre[ab_titre < eps] = eps
        ef_severe_uncond = 1 - (1 / (1 + np.exp(-k * (np.log10(ab_titre) - np.log10(ab_50_severe)))))
        ef_severe =  ef_severe_uncond / ef_infection

        return ef_severe

    def _update_ab_titre(self, vector_region, dr_vec, number_of_days):
        """Updates antibody titres"""

        if hasattr(vector_region, 'days_since_last_dose'):
            days_since_last = np.minimum(vector_region.days_since_last_infection,
                                         vector_region.days_since_last_dose)
        else:
            days_since_last = vector_region.days_since_last_infection
        indexes = np.argwhere(days_since_last < number_of_days + 1).flatten()
        valid_days_since_last = days_since_last[indexes]
        valid_days_since_last[valid_days_since_last > len(dr_vec) - 1] = len(dr_vec) - 1

        vector_region.ab_titre[indexes] =\
            vector_region.ab_titre[indexes] + dr_vec[valid_days_since_last]

    def _recovery(self, vector_region, new_recoveries, day):
        """Update antibody titres following an infection"""

        std10_infection = 0.44
        max_ab = 8
        mu_ab_inf = np.array([13 / 94, 223 / 94]) # may depend on day

        inf_number = vector_region.infection_number[new_recoveries]
        inf_number[inf_number > len(mu_ab_inf) - 1] = len(mu_ab_inf) - 1

        zdose = np.power(10, np.random.normal(np.log10(mu_ab_inf[inf_number]), std10_infection,
                                            new_recoveries.size))

        vector_region.ab_titre[new_recoveries] =\
            np.log(np.exp(vector_region.ab_titre[new_recoveries]) + zdose)
        vector_region.ab_titre[vector_region.ab_titre > max_ab] = max_ab
        vector_region.infection_number[new_recoveries] += 1
        vector_region.days_since_last_infection[new_recoveries] = 0
        vector_region.current_strain[new_recoveries] = -1

    def _load_contact_matrix(self, vector_region):
        """Loads contact matrix"""

        if (self.age_mixing_matrices_directory is not None) and self.sir_rescaling:

            # Discrete age
            number_of_age_groups = self.number_of_age_groups
            discrete_age = np.floor_divide(vector_region.age, self.age_group_interval)
            discrete_age[discrete_age > number_of_age_groups - 1] = number_of_age_groups - 1
            discrete_age = discrete_age.astype(np.int64)
            pop_by_disc_age = np.zeros(number_of_age_groups, dtype=np.int64)
            for a in range(number_of_age_groups):
                count = np.count_nonzero(vector_region.age == a)
                pop_by_disc_age[a] = count

            # Contact matrix
            filepath = self.age_mixing_matrices_directory + '/' + vector_region.other_name + '.csv'
            contact_matrix = np.genfromtxt(filepath, delimiter=',', dtype=np.float64)
            contact_matrix = np.vstack((contact_matrix, contact_matrix[15]))
            contact_matrix = np.hstack((contact_matrix, contact_matrix[:,15][:, np.newaxis] *\
                (pop_by_disc_age[16] / (pop_by_disc_age[15] + pop_by_disc_age[16]))))
            contact_matrix[:,15] = contact_matrix[:,15] *\
                (pop_by_disc_age[15] / (pop_by_disc_age[15] + pop_by_disc_age[16]))
            contact_matrix = np.multiply(contact_matrix, pop_by_disc_age[:, np.newaxis])
            contact_matrix = (contact_matrix + contact_matrix.T) / 2
            contact_matrix = np.divide(contact_matrix, pop_by_disc_age[:, np.newaxis])
            contact_matrix = contact_matrix / pop_by_disc_age[np.newaxis, :]

            vector_region.number_of_subpopulations = number_of_age_groups
            vector_region.subpopulation_index = discrete_age
            vector_region.subpopulation_mixing_matrix = contact_matrix

        else:

            vector_region.number_of_subpopulations = 1
            vector_region.subpopulation_index =\
                np.zeros((vector_region.number_of_agents), dtype=np.int64)
            vector_region.subpopulation_mixing_matrix = np.ones((1, 1), dtype=np.float64)
