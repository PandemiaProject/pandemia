"""This file creates a world."""

import logging
import csv
import numpy as np
import shapefile as shp
from scipy.spatial import distance
import random

from ..agent import Agent 
from ..location import Location
from ..region import Region
from .. import World
from . import WorldFactory

log = logging.getLogger('world_factory')

class SIRWorldFactory(WorldFactory):
    """Creates a simple world consisting of a single region, with one activity and one location.
    Using the sir_rescaling option in the default health model, agents in this world will mix
    homogeneously, as in the SIR model."""

    def __init__(self, config, clock, scale_factor):
        """Create agents and locations according to the population density map given"""

        self.config = config
        self.clock = clock
        self.scale_factor = scale_factor

        self.region_name = self.config['region_name']
        self.N           = self.config['N']

    def get_world(self) -> World:

        log.info("Creating world...")

        world = World(self.scale_factor)

        self._create_region(self.clock.ticks_in_week, world)

        number_of_regions = len(world.regions)

        assert number_of_regions == 1

        world.travel_matrix =\
            np.zeros((number_of_regions, number_of_regions), dtype=np.int64)

        return world

    def _create_region(self, ticks_in_week, world):
        """Initialize regions"""

        # Initialize regions using population data
        activities = ["Default"]
        world.regions = []
        new_location = Location(self.region_name, (0.0, 0.0))
        locations = [new_location]
        number_of_agents = max(int(self.N * self.scale_factor), 1)
        agents = []
        for _ in range(number_of_agents):
            age = random.randrange(100)
            new_agent = Agent(age)
            new_agent.weekly_routine = ["Default" for _ in range(ticks_in_week)]
            new_agent.activity_locations = {"Default": [new_location]}
            new_agent.activity_location_weights = {"Default": [1.0]}
            agents.append(new_agent)
        new_region = Region(0, self.region_name, activities, agents, locations)
        world.regions.append(new_region)
        world.number_of_regions = len(world.regions)

        log.info("Created world with %d agents", number_of_agents)
