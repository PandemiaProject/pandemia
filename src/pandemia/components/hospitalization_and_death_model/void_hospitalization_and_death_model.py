"""Void hospitalization and death model"""

import logging

from ..hospitalization_and_death_model import HospitalizationAndDeathModel

log = logging.getLogger("void_hospitalization_and_death_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class VoidHospitalizationAndDeathModel(HospitalizationAndDeathModel):
    """Void model of agent hospitalization and death"""

    def __init__(self, config):
        """Initial agent hospitalization and death"""
        super().__init__(config)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        pass

    def initial_conditions(self, vector_region):

        pass

    def dynamics(self, vector_region):
        """Changes agent hospitalization and death"""

        pass
