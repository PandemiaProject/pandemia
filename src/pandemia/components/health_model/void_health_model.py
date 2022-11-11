"""Void health model"""

import logging
import copy
import numpy as np

from .components.health_model import HealthModel

log = logging.getLogger("void_health_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class VoidHealthModel(HealthModel):
    """Void model of agent health"""

    def __init__(self, config, scale_factor, clock):
        """Initial agent health"""
        super().__init__(config, scale_factor)

        self.number_of_strains = 0
        self.immunity_length = 1
        self.immunity_period_ticks = 1

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        pass

    def initial_conditions(self, vector_region):

        pass

    def dynamics(self, vector_region, t):
        """Changes agent health"""

        pass
