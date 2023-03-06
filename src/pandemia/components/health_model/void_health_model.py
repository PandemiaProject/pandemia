"""Void health model"""

import logging
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

        self.number_of_strains = 1
        self.number_of_rho_immunity_outcomes = 1

        self.immunity_period_days = 365
        self.immunity_length = (clock.simulation_length_days // self.immunity_period_days) + 1
        self.immunity_period_ticks = self.immunity_period_days * clock.ticks_in_day
        self.facemask_transmission_multiplier = 1.0
        self.beta = np.array([0.0], dtype=np.float64)
        self.mutation_matrix = np.asarray([[1.0]], dtype=np.float64).flatten()

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        number_of_agents = vector_region.number_of_agents

        vector_region.current_disease = np.zeros((number_of_agents), dtype=np.float64)
        vector_region.current_strain =  np.full((number_of_agents), -1, dtype=np.int32)
        vector_region.current_infectiousness =  np.zeros((number_of_agents), dtype=np.float64)
        vector_region.requesting_immunity_update = np.zeros((number_of_agents), dtype=np.int32)
        vector_region.rho_immunity_failure_values = np.ones((number_of_agents,
                                                             self.number_of_strains,
                                                             self.immunity_length,
                                                             self.number_of_rho_immunity_outcomes),
                                                             dtype=np.float64)
        vector_region.sigma_immunity_failure_values = np.ones((number_of_agents,
                                                               self.number_of_strains,
                                                               self.immunity_length),
                                                               dtype=np.float64)
        vector_region.current_rho_immunity_failure = np.ones((number_of_agents,
                                                              self.number_of_strains,
                                                              self.number_of_rho_immunity_outcomes),
                                                              dtype=np.float64)
        vector_region.current_sigma_immunity_failure = np.ones((number_of_agents,
                                                                self.number_of_strains),
                                                                dtype=np.float64)
        vector_region.infection_event = np.full((number_of_agents), -1, dtype=np.int32)

    def initial_conditions(self, vector_region):

        pass

    def dynamics(self, vector_region, t):
        """Changes agent health"""

        pass

    def infect_wrapper(self, vector_region, t):
        """Changes to agent health"""

        pass

    def update(self, vector_region, t):
        """Updates health functions"""

        pass