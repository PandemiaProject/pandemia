"""Safir vaccination model"""

import logging

import numpy as np

from pandemia.components.vaccination_model import VaccinationModel

log = logging.getLogger("safir_vaccination_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class SafirVaccinationModel(VaccinationModel):
    """A version of the Safir3 vaccination model:

    To see the Safir3 model, follow the link: https://github.com/mrc-ide/safir3.

    Parameters:
    -----------
    config : Config
        A Pandemia Config object. A sub-config of the full config, containing the data used to
        configure this component.
    clock : Clock
        A Pandemia Clock object, discretizing the day.
    number_of_strains : int
        Not used in the Safir model.
    number_of_rho_immunity_outcomes : int
        Not used in the Safir model.
    immunity_length : int
        Not used in the Safir model. 
    immunity_period_ticks : int
        Not used in the Safir model.
    """

    def __init__(self, config, clock, number_of_strains, number_of_rho_immunity_outcomes,
                 immunity_length, immunity_period_ticks):
        """Initialize component"""
        super().__init__(config)

        self.clock = clock
        self.number_of_days = self.clock.simulation_length_days

        self.number_of_vaccines = 1
        self.age_groups = [0]

        # Vaccination strategy
        self.days_between_first_and_second_doses = 21
        self.number_of_age_groups = 17                                                              # TODO: must be same as in health model, enforce this...
        self.number_of_steps = 17
        coverage = 1.0
        self.strategy_matrix =\
            np.flip(np.triu(np.full((self.number_of_steps, self.number_of_age_groups),
                                     coverage, dtype=np.float64)), axis=0)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        number_of_agents = vector_region.number_of_agents

        vector_region.first_doses_total = 0
        vector_region.third_doses_total = 0

        vector_region.dose_number =\
            np.zeros(number_of_agents, dtype=np.int64)
        vector_region.days_since_last_dose =\
            np.full(number_of_agents, self.number_of_days + 1, dtype=np.int64)
        vector_region.day_for_second_dose =\
            np.full(number_of_agents, -1, dtype=np.int64)

    def initial_conditions(self, vector_region):
        """Initial vaccination"""

        # Vaccination order
        order = self._get_order(vector_region, self.strategy_matrix)
        vector_region.first_dose_order = order
        vector_region.third_dose_order = order

    def dynamics(self, vector_region, day, ticks_in_day):
        """Changes to agent health"""

        # Vaccinations (first doses)
        num_new_first_doses = 0
        low = min(vector_region.first_doses_total, vector_region.first_dose_order.size)
        high = min(vector_region.first_doses_total + num_new_first_doses,
                    vector_region.first_dose_order.size)
        new_doses = vector_region.first_dose_order[low:high]
        self._vaccinate(vector_region, new_doses, 0)
        vector_region.day_for_second_dose[new_doses] =\
            day + self.days_between_first_and_second_doses
        vector_region.first_doses_total += new_doses.size

        # Vaccinations (second doses)
        new_doses = np.argwhere(vector_region.day_for_second_dose == day).flatten()
        self._vaccinate(vector_region, new_doses, 1)

        # Vaccinations (third doses)
        num_new_third_doses = 0
        low = min(vector_region.third_doses_total, vector_region.third_dose_order.size)
        high = min(vector_region.third_doses_total + num_new_third_doses,
                    vector_region.third_dose_order.size)
        new_doses = vector_region.third_dose_order[low:high]
        self._vaccinate(vector_region, new_doses, 2)
        vector_region.third_doses_total += new_doses.size

        # Update days since last dose
        indexes = np.argwhere(vector_region.days_since_last_dose <
                                self.number_of_days + 1).flatten()
        vector_region.days_since_last_dose[indexes] += 1

    def _get_order(self, vector_region, strategy_matrix):
        """Determines the default order in which individuals are vaccinated"""

        age_groups =\
            np.zeros((self.number_of_age_groups, vector_region.number_of_agents), dtype=np.int64)
        if hasattr(vector_region, 'subpopulation_index'):
            for a in range(self.number_of_age_groups):
                age_groups[a][vector_region.subpopulation_index == a] = 1
        age_group_sizes = np.sum(age_groups, axis = 1)

        age_groups_copy = np.copy(age_groups)

        order = []
        for row in range(self.number_of_steps):
            sub_order = []
            for col in range(self.number_of_age_groups):
                coverage = strategy_matrix[row][col]
                age_group_size = age_group_sizes[col]
                num_to_vaccinate =\
                    max(int(np.sum(age_groups_copy[col]) - (age_group_size * (1 - coverage))), 0)
                if num_to_vaccinate > 0:
                    age_group_indexes = np.argwhere(age_groups_copy[col] == 1).flatten().tolist()
                    indexes_to_vaccinate =\
                        vector_region.prng.random_sample(age_group_indexes, num_to_vaccinate)
                    sub_order += indexes_to_vaccinate
                    age_groups_copy[col][indexes_to_vaccinate] = 0
            vector_region.prng.random_shuffle(sub_order)
            order += sub_order

        order = np.array(order, dtype=np.int64)

        return order

    def _vaccinate(self, vector_region, new_doses, dose):
        """Update antibody titres following a dose of vaccine"""

        std10 = 0.44
        max_ab = 8

        mu_ab_list = {
            "Pfizer": [13 / 94, 223 / 94],
            "AstraZeneca": [1 / 59, 32 / 59],
            "Sinovac": [28 / 164, 28 / 164],
            "Moderna": [((185 + 273) / 2) / 321, 654 / 158]
        }

        vaccine = 'Pfizer'

        mu_ab = np.array(mu_ab_list[vaccine])
        dose = min(dose, mu_ab.size - 1)

        zdose = np.power(10, np.random.normal(np.log10(mu_ab[dose]), std10, new_doses.size))

        vector_region.ab_titre[new_doses] =\
            np.log(np.exp(vector_region.ab_titre[new_doses]) + zdose)
        vector_region.ab_titre[vector_region.ab_titre > max_ab] = max_ab

        vector_region.dose_number[new_doses] += 1
        vector_region.days_since_last_dose[new_doses] = 0
