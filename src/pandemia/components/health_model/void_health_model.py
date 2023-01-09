"""Void health model"""

import logging
import copy
import numpy as np

from ..health_model import HealthModel

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
        self.mutation_matrix = [[1.0]]
        
        # Copied from Default model
        self.facemask_transmission_multiplier = 1.0


    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        # This attributes are currently initialised in DefaultHealthModel and
        # copied here
        vector_region.current_disease = np.zeros((vector_region.number_of_agents), dtype=np.float64)
        vector_region.current_strain = np.full((vector_region.number_of_agents), -1, dtype=np.int64)

    def initial_conditions(self, vector_region):

        pass

    def dynamics(self, vector_region, t):
        """Changes agent health"""

        pass


    def update(self, vector_region, t):
        """Updates health functions"""

        pass