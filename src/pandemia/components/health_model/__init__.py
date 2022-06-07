"""Manages health model"""

from pandemia.component import Component

class HealthModel(Component):
    """Defines the health of agents"""

    def __init__(self, config, scale_factor):

        super().__init__(config)
