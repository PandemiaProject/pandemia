"""Manages seasonal effects"""

from ...component import Component

class SeasonalEffectsModel(Component):
    """Defines the seasonal effects model"""

    def __init__(self, config):

        super().__init__(config)
