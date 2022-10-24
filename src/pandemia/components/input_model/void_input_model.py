"""Void input model, specifiying a model with no policy interventions"""

import logging

from pandemia.components.input_model import InputModel

log = logging.getLogger("void_input_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class VoidInputModel(InputModel):
    """Default model of input, specifiying a model with no policy interventions"""

    def __init__(self, config, scale_factor, clock, number_of_regions,
                 number_of_vaccines, age_groups):
        """Initialize component"""
        super().__init__(config, scale_factor)

    def new_input(self, policy):
        """Set new input"""

        pass

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        pass

    def initial_conditions(self, vector_region):

        pass

    def dynamics(self, vector_region, day):
        """Changes to input"""

        pass
