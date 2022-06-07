"""Manages hospitalization and death"""

from pandemia.component import Component

class HospitalizationAndDeathModel(Component):
    """Defines the hospitalization and death model"""

    def __init__(self, config):

        super().__init__(config)
