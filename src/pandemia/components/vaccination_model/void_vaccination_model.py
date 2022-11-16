"""Void vaccination model"""

import logging

from ..vaccination_model import VaccinationModel

log = logging.getLogger("void_vaccination_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class VoidVaccinationModel(VaccinationModel):
    """Void model of vaccination"""

    def __init__(self, config, clock, number_of_strains, immunity_length, immunity_period_ticks):
        """Initial vaccination"""
        super().__init__(config)

        self.number_of_vaccines = 0
        self.age_groups = [0]

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        pass

    def initial_conditions(self, vector_region):
        """Initial vaccination"""

        pass

    def dynamics(self, vector_region, day, ticks_in_day):
        """Changes vaccination"""

        pass
