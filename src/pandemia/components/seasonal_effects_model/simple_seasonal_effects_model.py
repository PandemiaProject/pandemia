"""Simple seasonal effects model"""

import logging
import numpy as np
from datetime import timedelta

from pandemia.components.seasonal_effects_model import SeasonalEffectsModel

log = logging.getLogger("simple_seasonal_effects_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class SimpleSeasonalEffectsModel(SeasonalEffectsModel):
    """Simple model of seasonal effects"""

    def __init__(self, config, clock, enable_ctypes):
        """Initialize component"""
        super().__init__(config)

        self.epoch = clock.epoch
        self.simulation_length_days = clock.simulation_length_days

        self.seasonal_multiplier_by_region = config['seasonal_multiplier_by_region_by_month']

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        vector_region.region_transmission_multiplier =\
            np.ones((self.simulation_length_days), dtype=float)
        vector_region.current_region_transmission_multiplier = 1.0

    def initial_conditions(self, vector_region):
        """Initial seasonal effect"""

        for day in range(self.simulation_length_days):
            date = self.epoch + timedelta(days=day)
            vector_region.region_transmission_multiplier[day] =\
                float(self.seasonal_multiplier_by_region[vector_region.id][date.month - 1])

        vector_region.current_region_transmission_multiplier =\
            vector_region.region_transmission_multiplier[0]

    def dynamics(self, vector_region, day):
        """Changes to seasonal effects"""

        vector_region.current_region_transmission_multiplier =\
            vector_region.region_transmission_multiplier[day]
