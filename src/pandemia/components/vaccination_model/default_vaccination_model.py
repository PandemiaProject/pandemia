"""Default vaccination model"""

import logging

import numpy as np

from ctypes import c_void_p, c_int, cdll

from ..vaccination_model import VaccinationModel
log = logging.getLogger("default_vaccination_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultVaccinationModel(VaccinationModel):
    """Default model of vaccination. This supports multiple vaccines, acting on multiple strains of
    the pathogen. Vaccines update individual immunity according to the same mechanism used to update
    immunity in the default health model."""

    def __init__(self, config, clock, number_of_strains, number_of_rho_immunity_outcomes,
                 immunity_length, immunity_period_ticks):
        """Initialize component"""
        super().__init__(config)

        self.dynamics_vaccination = self.lib.dynamics_vaccination

        self.dynamics_vaccination.restype = None

        self.booster_waiting_time_days = config['booster_waiting_time_days']
        self.age_groups                = config['age_groups']
        self.number_of_age_groups      = len(self.age_groups)
        self.vaccine_hesitancy         = config['vaccine_hesitancy']
        self.vaccines_config           = config['vaccines']

        self.number_of_strains         = number_of_strains
        self.immunity_length           = immunity_length
        self.immunity_period_ticks     = immunity_period_ticks

        self.vaccine_names             = list(self.vaccines_config.keys())
        self.number_of_vaccines        = len(self.vaccines_config)

        self.number_of_strains = self._get_number_of_strains()
        assert self.number_of_strains == number_of_strains,\
            "Number of strains in vaccination model must match health model"

        self.number_of_rho_immunity_outcomes = self._get_number_of_rho_immunity_outcomes()
        assert self.number_of_rho_immunity_outcomes == number_of_rho_immunity_outcomes,\
            "Number of rho immunity outcomes in vaccination model must match health model"

        self.max_vaccine_length_immunity     = self._get_max_vaccine_length_immunity()

        self.vaccine_rho_immunity_failure_values =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains,
                      self.max_vaccine_length_immunity,
                      self.number_of_rho_immunity_outcomes), dtype=np.float64)

        self.vaccine_rho_immunity_failure_partitions =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains,
                      self.max_vaccine_length_immunity), dtype=np.int64)

        self.vaccine_rho_immunity_failure_lengths =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains), dtype=np.int64)

        self.vaccine_sigma_immunity_failure_values =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains,
                      self.max_vaccine_length_immunity), dtype=np.float64)

        self.vaccine_sigma_immunity_failure_partitions =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains,
                      self.max_vaccine_length_immunity), dtype=np.int64)

        self.vaccine_sigma_immunity_failure_lengths =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains), dtype=np.int64)

        # Create vaccines
        self._get_vaccines(clock.ticks_in_day)

        # Flatten arrays
        self.vaccine_rho_immunity_failure_values =\
            self.vaccine_rho_immunity_failure_values.flatten()
        self.vaccine_rho_immunity_failure_partitions =\
            self.vaccine_rho_immunity_failure_partitions.flatten()
        self.vaccine_rho_immunity_failure_lengths =\
            self.vaccine_rho_immunity_failure_lengths.flatten()
        self.vaccine_sigma_immunity_failure_values =\
            self.vaccine_sigma_immunity_failure_values.flatten()
        self.vaccine_sigma_immunity_failure_partitions =\
            self.vaccine_sigma_immunity_failure_partitions.flatten()
        self.vaccine_sigma_immunity_failure_lengths =\
            self.vaccine_sigma_immunity_failure_lengths.flatten()

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        number_of_agents = vector_region.number_of_agents
        number_of_vaccines = self.number_of_vaccines

        vector_region.number_of_vaccination_age_groups = 0
        vector_region.vaccination_age_group = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.num_to_vaccinate =\
            np.zeros((vector_region.number_of_vaccination_age_groups,
                      number_of_vaccines), dtype=np.int64)
        vector_region.most_recent_first_dose = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.vaccine_hesitant = np.zeros((number_of_agents), dtype=np.int64)

    def initial_conditions(self, vector_region):
        """Initial vaccination"""

        # Determine who belongs to each vaccination age group and set other initial data
        for n in range(vector_region.number_of_agents):
            age = vector_region.age[n]
            age_group_index = self._determine_age_group_index(age, self.age_groups)
            vector_region.vaccination_age_group[n] = age_group_index
            vector_region.most_recent_first_dose[n] = - self.booster_waiting_time_days
            if vector_region.prng.random_float(1.0) < self.vaccine_hesitancy[age_group_index]:
                vector_region.vaccine_hesitant[n] = 1

    def dynamics(self, vector_region, day, ticks_in_day):
        """Changes to agent health"""

        if np.sum(vector_region.num_to_vaccinate) > 0:
            vector_region.num_to_vaccinate = vector_region.num_to_vaccinate.flatten()
            self.dynamics_vaccination(
                c_int(day),
                c_int(ticks_in_day),
                c_int(vector_region.number_of_agents),
                c_int(self.number_of_strains),
                c_int(self.immunity_length),
                c_int(self.number_of_rho_immunity_outcomes),
                c_int(self.number_of_age_groups),
                c_int(self.number_of_vaccines),
                c_int(self.max_vaccine_length_immunity),
                c_int(self.booster_waiting_time_days),
                c_int(self.immunity_period_ticks),
                c_int(vector_region.id),
                c_void_p(self.vaccine_rho_immunity_failure_values.ctypes.data),
                c_void_p(self.vaccine_rho_immunity_failure_partitions.ctypes.data),
                c_void_p(self.vaccine_rho_immunity_failure_lengths.ctypes.data),
                c_void_p(self.vaccine_sigma_immunity_failure_values.ctypes.data),
                c_void_p(self.vaccine_sigma_immunity_failure_partitions.ctypes.data),
                c_void_p(self.vaccine_sigma_immunity_failure_lengths.ctypes.data),
                c_void_p(vector_region.num_to_vaccinate.ctypes.data),
                c_void_p(vector_region.rho_immunity_failure_values.ctypes.data),
                c_void_p(vector_region.sigma_immunity_failure_values.ctypes.data),
                c_void_p(vector_region.current_region.ctypes.data),
                c_void_p(vector_region.current_disease.ctypes.data),
                c_void_p(vector_region.current_strain.ctypes.data),
                c_void_p(vector_region.requesting_immunity_update.ctypes.data),
                c_void_p(vector_region.most_recent_first_dose.ctypes.data),
                c_void_p(vector_region.vaccine_hesitant.ctypes.data),
                c_void_p(vector_region.random_state.ctypes.data)
            )

    def _get_number_of_strains(self):
        """Determines number of strains in config"""

        num_strains = len(self.vaccines_config[self.vaccine_names[0]]['rho_immunity_failure'])

        for vaccine_name in self.vaccine_names:
            assert len(self.vaccines_config[vaccine_name]['rho_immunity_failure']) == num_strains
            assert len(self.vaccines_config[vaccine_name]['sigma_immunity_failure']) == num_strains

        return num_strains

    def _get_number_of_rho_immunity_outcomes(self):
        """Determines number of rho immunity outcomes in config"""

        num_rho_outcomes =\
            len(self.vaccines_config[self.vaccine_names[0]]['rho_immunity_failure'][0][1][0])

        for vaccine_name in self.vaccine_names:
            for s in range(self.number_of_strains):
                items = len(self.vaccines_config[vaccine_name]['rho_immunity_failure'][s][1])
                for i in range(items):
                    assert len(self.vaccines_config[vaccine_name]['rho_immunity_failure'][s][1][i])\
                           == num_rho_outcomes

        return num_rho_outcomes

    def _get_max_vaccine_length_immunity(self):
        """Determine maximum length of vaccine immunity outcomes in config"""

        max_vaccine_length_immunity = 0
        for vaccine_name in self.vaccines_config:
            for s in range(self.number_of_strains):
                for immunity_failure_type in ['rho_immunity_failure', 'sigma_immunity_failure']:
                    partition_data = self.vaccines_config[vaccine_name][immunity_failure_type][s][0]
                    values_data = self.vaccines_config[vaccine_name][immunity_failure_type][s][1]
                    assert len(partition_data) == len(values_data)
                    if len(partition_data) > max_vaccine_length_immunity:
                        max_vaccine_length_immunity = len(partition_data)

        return max_vaccine_length_immunity

    def _get_vaccines(self, ticks_in_day):
        """Create vaccines"""

        # Build vaccines
        for vaccine_name in self.vaccines_config:
            id = self.vaccine_names.index(vaccine_name)
            f = 'rho_immunity_failure'
            for s in range(self.number_of_strains):
                partition_data = self.vaccines_config[vaccine_name][f][s][0]
                values_data = self.vaccines_config[vaccine_name][f][s][1]
                partition = [int(part * ticks_in_day) for part in partition_data]
                values = [[float(val1) for val1 in val2] for val2 in values_data]
                length = len(partition)
                self.vaccine_rho_immunity_failure_lengths[id][s] = length
                for i in range(length):
                    self.vaccine_rho_immunity_failure_partitions[id][s][i] = partition[i]
                    for r in range(self.number_of_rho_immunity_outcomes):
                        self.vaccine_rho_immunity_failure_values[id][s][i][r] = values[i][r]
                for i in range(length, self.max_vaccine_length_immunity):
                    self.vaccine_rho_immunity_failure_partitions[id][s][i] = -1
                    for r in range(self.number_of_rho_immunity_outcomes):
                        self.vaccine_rho_immunity_failure_values[id][s][i][r] = 1.0
            f = 'sigma_immunity_failure'
            for s in range(self.number_of_strains):
                partition_data = self.vaccines_config[vaccine_name][f][s][0]
                values_data = self.vaccines_config[vaccine_name][f][s][1]
                partition = [int(part * ticks_in_day) for part in partition_data]
                values = [float(val) for val in values_data]
                length = len(partition)
                self.vaccine_sigma_immunity_failure_lengths[id][s] = length
                for i in range(length):
                    self.vaccine_sigma_immunity_failure_partitions[id][s][i] = partition[i]
                    self.vaccine_sigma_immunity_failure_values[id][s][i] = values[i]
                for i in range(length, self.max_vaccine_length_immunity):
                    self.vaccine_sigma_immunity_failure_partitions[id][s][i] = -1
                    self.vaccine_sigma_immunity_failure_values[id][s][i] = 1.0

    def _determine_age_group_index(self, age, age_groups):
        """Determine to which age group the agent belongs"""

        if age >= age_groups[-1]:
            return len(age_groups) - 1
        elif age < age_groups[0]:
            return 0
        else:
            i = 0
            while age < age_groups[i] or age >= age_groups[i+1]:
                i = i + 1
            return i
