"""Void testing and contact tracing model"""

import logging

from .components.testing_and_contact_tracing_model import TestingAndContactTracingModel

log = logging.getLogger("void_testing_and_contact tracing_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class VoidTestingAndContactTracingModel(TestingAndContactTracingModel):
    """Void model of testing and contact tracing"""

    def __init__(self, config):
        """Initial testing and contact tracing"""
        super().__init__(config)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        pass

    def initial_conditions(self, vector_region):

        pass

    def dynamics(self, vector_region, day):
        """Changes testing and contact tracing"""

        pass
