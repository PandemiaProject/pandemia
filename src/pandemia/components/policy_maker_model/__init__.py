"""Manages policy interventions in the model"""

from ...component import Component

class PolicyMakerModel(Component):
    """Defines the policy maker in the model"""

    def __init__(self, config, scale_factor):

        super().__init__(config)

    def new_policy(self):
        """Set new policy"""
        pass
