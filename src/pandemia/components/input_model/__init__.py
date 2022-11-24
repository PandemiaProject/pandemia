"""Manages input to the model"""

from ...component import Component

class InputModel(Component):
    """Defines the input to the model"""

    def __init__(self, config, scale_factor):

        super().__init__(config)

    def new_input(self):
        """Set new input"""
        pass
