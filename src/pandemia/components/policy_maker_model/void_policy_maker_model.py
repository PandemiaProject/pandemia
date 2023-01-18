"""Void policy model, specifiying a model with no policy interventions"""

import logging

from ..policy_maker_model import PolicyMakerModel

log = logging.getLogger("void_policy_maker_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class VoidPolicyMakerModel(PolicyMakerModel):
    """Default model of a policy maker, specifiying a model with no policy interventions"""

    def __init__(self, config, scale_factor, clock, number_of_regions,
                 number_of_vaccines, age_groups):
        """Initialize component"""
        super().__init__(config, scale_factor)

    def new_policy(self, policy):
        """Set new policy"""

        pass

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        pass

    def initial_conditions(self, vector_region):

        pass

    def dynamics(self, vector_region, day):
        """Changes to policy"""

        pass
