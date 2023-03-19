"""Void regional mixing model"""

import logging
import numpy as np

from ..travel_model import TravelModel

log = logging.getLogger("void_travel_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class VoidTravelModel(TravelModel):
    """Void model of regional mixing"""

    def __init__(self, config, scale_factor, number_of_strains,
                 number_of_regions):
        """Initial regional mixing"""
        super().__init__(config, scale_factor)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        vector_region.current_region = np.full((vector_region.number_of_agents),
                                                vector_region.id, dtype=np.int64)

    def initial_conditions(self, sim):

        pass

    def dynamics(self, sim, day, ticks_in_day,
                 facemask_transmission_multiplier,
                 mutation_matrix, enable_parallel,
                 num_jobs, vector_region_batches):
        """Changes regional mixing"""

        pass
