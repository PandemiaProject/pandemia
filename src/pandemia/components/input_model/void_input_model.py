"""Void input model"""

import logging

from pandemia.components.input_model import InputModel

log = logging.getLogger("void_input_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class VoidInputModel(InputModel):
    """Simple model of input"""

    def __init__(self, config, scale_factor, simulation_length_days, number_of_regions,
                 number_of_vaccines, age_groups, enable_ctypes):
        """Initialize component"""
        super().__init__(config, scale_factor)

    def new_input(self, input_arrays):
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
