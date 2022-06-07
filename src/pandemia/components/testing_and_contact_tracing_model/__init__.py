"""Manages testing and contact tracing"""

from pandemia.component import Component

class TestingAndContactTracingModel(Component):
    """Defines the testing and contact tracing model"""

    def __init__(self, config):

        super().__init__(config)
