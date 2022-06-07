"""Manages regional mixing"""

from pandemia.component import Component

class RegionalMixingModel(Component):
    """Defines the regional mixing model"""

    def __init__(self, config, scale_factor):

        super().__init__(config)
