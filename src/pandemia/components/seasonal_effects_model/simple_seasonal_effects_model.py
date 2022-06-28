"""Simple seasonal effects model"""

from collections import defaultdict
import logging
import numpy as np
import csv
from datetime import timedelta

from pandemia.components.seasonal_effects_model import SeasonalEffectsModel

log = logging.getLogger("simple_seasonal_effects_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class SimpleSeasonalEffectsModel(SeasonalEffectsModel):
    """Simple model of seasonal effects"""

    def __init__(self, config, vector_world, clock, enable_ctypes):
        """Initialize component"""
        super().__init__(config)

        self.epoch = clock.epoch
        self.simulation_length_days = clock.simulation_length_days

        seasonality_config = config['seasonal_multiplier_by_region_by_month']

        self.seasonal_multiplier_by_region = defaultdict(list)

        if seasonality_config is None:
            for vector_region in vector_world.vector_regions:
                self.seasonal_multiplier_by_region[vector_region.id] = [1 for _ in range(12)]
        else:
            seasonal_multiplier_records = defaultdict(list)
            with open(seasonality_config, newline='') as csvfile:
                next(csvfile)
                region_data_seasonality = csv.reader(csvfile, delimiter=',')
                for row in region_data_seasonality:
                    iso = str(row[2])
                    if iso == 'GB':
                        iso = 'UK'
                    seasonal_multipliers = [float(row[3 + r]) for r in range(12)]
                    seasonal_multiplier_records[iso] = seasonal_multipliers
            names_to_ids = {r.name: r.id for r in vector_world.vector_regions}
            for vector_region in vector_world.vector_regions:
                if vector_region.name in seasonal_multiplier_records:
                    self.seasonal_multiplier_by_region[names_to_ids[vector_region.name]] =\
                        seasonal_multiplier_records[vector_region.name]

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
