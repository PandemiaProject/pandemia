"""Manages vaccination"""

from pandemia.component import Component

class VaccinationModel(Component):
    """Defines the vaccination model"""

    def __init__(self, config):

        super().__init__(config)
