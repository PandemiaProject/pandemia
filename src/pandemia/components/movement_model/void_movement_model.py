"""Void location selection."""

import logging
import numpy as np

from ..movement_model import MovementModel

log = logging.getLogger("void_movement_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class VoidMovementModel(MovementModel):
    """Void movement model"""

    def __init__(self, config):
        """Initial agent locations"""
        super().__init__(config)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        number_of_agents = vector_region.number_of_agents
        number_of_activities = vector_region.number_of_activities

        vector_region.requested_location_update = np.zeros((number_of_agents), dtype=np.int32)
        vector_region.requesting_location_update = np.zeros((number_of_agents), dtype=np.int32)
        vector_region.home_location = np.zeros((number_of_agents), dtype=np.int32)
        vector_region.current_quarantine = np.zeros((number_of_agents), dtype=np.int32)
        vector_region.current_location = np.zeros((number_of_agents), dtype=np.int32)
        vector_region.current_facemask = np.zeros((number_of_agents), dtype=np.int32)
        vector_region.requested_facemask_update = np.zeros((number_of_agents), dtype=np.int32)
        vector_region.requesting_facemask_update = np.zeros((number_of_agents), dtype=np.int32)
        vector_region.wears_facemask =\
            np.zeros((number_of_agents, number_of_activities), dtype=np.int32)

    def initial_conditions(self, vector_region, offset):

        pass

    def update(self, vector_region, t):
        """Updates movement functions"""

        pass

    def dynamics(self, vector_region, t, t_0, ticks_in_week):
        """Changes agent locations"""

        pass
