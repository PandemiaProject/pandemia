"""Manages regional mixing"""

from ...component import Component

class TravelModel(Component):
    """Defines the regional mixing model"""

    def __init__(self, config, scale_factor):

        super().__init__(config)
