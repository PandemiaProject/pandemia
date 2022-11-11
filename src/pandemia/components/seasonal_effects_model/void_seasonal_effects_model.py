"""Void seasonal effects model"""

import logging

from .components.seasonal_effects_model import SeasonalEffectsModel

log = logging.getLogger("void_seasonal_effects_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class VoidSeasonalEffectsModel(SeasonalEffectsModel):
    """Void model of seasonal effects"""

    def __init__(self, config, vector_world, clock):
        """Initial seasonal effects"""
        super().__init__(config)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        vector_region.current_region_transmission_multiplier = 1.0

    def initial_conditions(self, vector_region):

        pass

    def dynamics(self, vector_region, day):
        """Changes seasonal effects"""

        pass
