"""Default seasonal effects model"""

import logging
import numpy as np
import csv
from datetime import timedelta

from ..seasonal_effects_model import SeasonalEffectsModel

log = logging.getLogger("default_seasonal_effects_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultSeasonalEffectsModel(SeasonalEffectsModel):
    """Default model of seasonal effects. Each region is assigned a region_transmission_multiplier.
    Seasonal effects can implemented by varying this multiplier, which in this default model occurs
    for each region each month."""

    def __init__(self, config, vector_world, clock):
        """Initialize component"""
        super().__init__(config)

        self.epoch = clock.epoch
        self.simulation_length_days = clock.simulation_length_days

        seasonality_config            = config['seasonal_multiplier_by_region_by_month']
        self.out_of_season_multiplier = config['out_of_season_multiplier']

        self.seasonal_multiplier_by_region = {}

        if seasonality_config is None:
            for vector_region in vector_world.vector_regions:
                self.seasonal_multiplier_by_region[vector_region.id] = np.ones((12, ), dtype=np.float64)
        else:
            seasonal_multiplier_records = {}
            with open(seasonality_config, newline='') as csvfile:
                next(csvfile)
                region_data_seasonality = csv.reader(csvfile, delimiter=',')
                for row in region_data_seasonality:
                    iso = str(row[1])
                    seasonal_multipliers = []
                    for r in range(12):
                        if int(row[3 + r]) == 1:
                            multiplier = 1
                        else:
                            multiplier = 0
                        seasonal_multipliers.append(multiplier)
                    seasonal_multiplier_records[iso] = np.array(seasonal_multipliers, dtype=np.float64)
            names_to_ids = {r.name: r.id for r in vector_world.vector_regions}
            for vector_region in vector_world.vector_regions:
                if vector_region.name in seasonal_multiplier_records:
                    self.seasonal_multiplier_by_region[names_to_ids[vector_region.name]] =\
                        seasonal_multiplier_records[vector_region.name]

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        vector_region.region_transmission_multiplier =\
            np.ones((self.simulation_length_days), dtype=np.float64)
        vector_region.current_region_transmission_multiplier = 1.0

    def initial_conditions(self, vector_region):
        """Initial seasonal effect"""

        multipliers = self.seasonal_multiplier_by_region[vector_region.id]
        self.seasonal_multiplier_by_region[vector_region.id][multipliers == 0] =\
            self.out_of_season_multiplier

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
