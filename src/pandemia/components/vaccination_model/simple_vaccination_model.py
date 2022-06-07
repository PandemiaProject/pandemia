"""Simple vaccination model"""

import logging

import numpy as np

from ctypes import c_void_p, pointer, c_int, cdll

from pandemia.components.vaccination_model import VaccinationModel

log = logging.getLogger("simple_vaccination_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class SimpleVaccinationModel(VaccinationModel):
    """Simple model of vaccination"""

    def __init__(self, config, clock, number_of_strains, immunity_length, immunity_period_ticks,
                 enable_ctypes):
        """Initialize component"""
        super().__init__(config)

        self.enable_ctypes = enable_ctypes

        if self.enable_ctypes:

            lib = cdll.LoadLibrary("./src/pandemia/components/vaccination_model/"
                                   "simple_vaccination_model_functions.dll")

            self.dynamics_vaccination = lib.dynamics_vaccination

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

        self.number_of_rho_immunity_outcomes = self._get_number_of_rho_immunity_outcomes()
        self.max_vaccine_length_immunity     = self._get_max_vaccine_length_immunity()

        self.vaccine_rho_immunity_failure_values =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains,
                      self.max_vaccine_length_immunity,
                      self.number_of_rho_immunity_outcomes), dtype=float)

        self.vaccine_rho_immunity_failure_partitions =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains,
                      self.max_vaccine_length_immunity), dtype=int)

        self.vaccine_rho_immunity_failure_lengths =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains), dtype=int)

        self.vaccine_sigma_immunity_failure_values =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains,
                      self.max_vaccine_length_immunity), dtype=float)

        self.vaccine_sigma_immunity_failure_partitions =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains,
                      self.max_vaccine_length_immunity), dtype=int)

        self.vaccine_sigma_immunity_failure_lengths =\
            np.zeros((self.number_of_vaccines,
                      self.number_of_strains), dtype=int)

        # Create vaccines
        self._get_vaccines(clock.ticks_in_day)

        # If using ctypes, flatten arrays accordingly
        if self.enable_ctypes:
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
        vector_region.vaccination_age_group = np.zeros((number_of_agents), dtype=int)
        vector_region.num_to_vaccinate =\
            np.zeros((vector_region.number_of_vaccination_age_groups,
                      number_of_vaccines), dtype=int)
        vector_region.most_recent_first_dose = np.zeros((number_of_agents), dtype=int)
        vector_region.vaccine_hesitant = np.zeros((number_of_agents), dtype=int)

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
            if self.enable_ctypes:
                self.dynamics_c(vector_region, day, ticks_in_day)
            else:
                self.dynamics_python(vector_region, day, ticks_in_day)

    def dynamics_c(self, vector_region, day, ticks_in_day):
        """Change in vaccination"""

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
            c_void_p(vector_region.random_state.ctypes.data),
            pointer(vector_region.random_p)
        )

    def dynamics_python(self, vector_region, day, ticks_in_day):
        """Change in vaccination"""

        if np.sum(vector_region.num_to_vaccinate) == 0:
            return

        # Determine who is eligible to be vaccinated today
        eligible = {a : [] for a in range(self.number_of_age_groups)}
        for n in range(vector_region.number_of_agents):
            if vector_region.current_region[n] == vector_region.id\
                and vector_region.current_disease[n] < 1.0:
                if day - vector_region.most_recent_first_dose[n] >= self.booster_waiting_time_days:
                    if vector_region.current_strain[n] == -1:
                        if vector_region.vaccine_hesitant[n] == 0:
                            age_group = vector_region.vaccination_age_group[n]
                            eligible[age_group].append(n)

        # Vaccinate a sample of the eligible population and update their immunity functions
        t = day*ticks_in_day
        for age_group_index in range(self.number_of_age_groups):
            sample_size = min(len(eligible[age_group_index]),
                              sum(vector_region.num_to_vaccinate[age_group_index]))
            sample = vector_region.prng.random_sample(eligible[age_group_index], sample_size)
            v = 0
            num_vaccinated = 0
            for n in sample:
                if v < self.number_of_vaccines:
                    for s in range(self.number_of_strains):
                        # Determine new rho immunity
                        values_by_index = vector_region.rho_immunity_failure_values[n][s]
                        part = self._shift(self.vaccine_rho_immunity_failure_partitions[v][s], t)
                        values = self.vaccine_rho_immunity_failure_values[v][s]
                        length = self.vaccine_rho_immunity_failure_lengths[v][s]
                        new_values =\
                            self._rho_multiply_and_subsample(values_by_index, part, values, length)
                        for i in range(self.immunity_length):
                            for r in range(self.number_of_rho_immunity_outcomes):
                                vector_region.rho_immunity_failure_values[n][s][i][r] =\
                                    new_values[i][r]
                        # Determine new sigma immunity
                        values_by_index = vector_region.sigma_immunity_failure_values[n][s]
                        part = self._shift(self.vaccine_sigma_immunity_failure_partitions[v][s], t)
                        values = self.vaccine_sigma_immunity_failure_values[v][s]
                        length = self.vaccine_sigma_immunity_failure_lengths[v][s]
                        new_values =\
                            self._multiply_and_subsample(values_by_index, part, values, length)
                        for i in range(self.immunity_length):
                            vector_region.sigma_immunity_failure_values[n][s][i] =\
                                new_values[i]
                    vector_region.most_recent_first_dose[n] = day
                    num_vaccinated += 1
                if num_vaccinated >= vector_region.num_to_vaccinate[age_group_index][v]:
                    v += 1
                    num_vaccinated = 0

    def _shift(self, partition, t: int):
        """Shifts a partition by t"""

        shifted_partition = [p + t for p in partition]
        shifted_partition[0] = -1

        return shifted_partition

    def _evaluate(self, partition, values, length, t: int):
        """Evaluates a given step function"""

        if t < partition[1]:
            return values[0]
        else:
            i = length - 1
            while t < partition[i]:
                i = i - 1
            return values[i]

    def _multiply_and_subsample(self, values_by_index, partition, values, length):
        """Multiplies sigma immunity array by a step function"""

        new_values = []

        for index in range(self.immunity_length):
            t = (index + 1) * self.immunity_period_ticks
            new_value = values_by_index[index] * self._evaluate(partition, values, length, t)
            new_values.append(new_value)

        return new_values

    def _rho_evaluate(self, partition, values, length, rho_outcome, t: int):
        """Evaluates a given vector-valued step function"""

        if t < partition[1]:
            return values[0][rho_outcome]
        else:
            i = length - 1
            while t < partition[i]:
                i = i - 1
            return values[i][rho_outcome]

    def _rho_multiply_and_subsample(self, values_by_index, partition, values, length):
        """Multiplies rho immunity array by a vector-valued step function"""

        new_values = []

        for index in range(self.immunity_length):
            t = (index + 1) * self.immunity_period_ticks
            new_values_rho_outcome = []
            for r in range(self.number_of_rho_immunity_outcomes):
                new_value = values_by_index[index][r] *\
                            self._rho_evaluate(partition, values, length, r, t)
                new_values_rho_outcome.append(new_value)
            new_values.append(new_values_rho_outcome)

        return new_values

    def _get_number_of_rho_immunity_outcomes(self):
        """Determines number of rho immunity outcomes in config"""

        num_rho_outcomes = len(self.vaccines_config[self.vaccine_names[0]])
        for vaccine_name in self.vaccines_config:
            assert len(self.vaccines_config[vaccine_name]) == num_rho_outcomes

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
