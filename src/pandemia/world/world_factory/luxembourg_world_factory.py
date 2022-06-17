"""This file creates the world, adding the map and a population of agents."""

import logging
import csv
import numpy as np
import shapefile as shp
from scipy.spatial import distance

from pandemia.world.agent import Agent
from pandemia.world.location import Location
from pandemia.world.region import Region
from pandemia.world import World
from pandemia.world.world_factory import WorldFactory

log = logging.getLogger('world_factory')

class LuxembourgWorldFactory(WorldFactory):
    """Reads a DensityMap and generates a world based on the densities indicated therein."""

    def __init__(self, config, clock, scale_factor):
        """Create agents and locations according to the population density map given"""

        self.config = config
        self.clock = clock
        self.scale_factor = scale_factor

        self.luxembourg_data_path = self.config['luxembourg_data_path']

    def get_world(self) -> World:

        log.info("Creating world...")

        world = World(self.scale_factor)

        self._create_luxembourg(world)

        number_of_regions = len(world.regions)

        assert number_of_regions == 1

        world.travel_matrix =\
            np.zeros((number_of_regions, number_of_regions), dtype=int)
        world.contacts_matrix =\
            np.full((number_of_regions, number_of_regions),
                    self.contacts, dtype=int)
        np.fill_diagonal(world.contacts_matrix, 0)
        contact_hours_matrix =\
            np.full((number_of_regions, number_of_regions),
                    self.contact_hours, dtype=int)
        np.fill_diagonal(contact_hours_matrix, 0)
        world.contact_ticks_matrix =\
            (self.clock.ticks_in_hour * contact_hours_matrix).astype(int)

        return world

    def _create_luxembourg(self, world):
        """Creates world for example one"""

        ticks_in_week = self.clock.ticks_in_week
        regions_data_path = self.regions_data_path

        # Initialize regions using population data
        activities = ["Default"]
        world.regions = []
        total_number_of_agents = 0
        with open(regions_data_path, newline='') as csvfile:
            region_data = csv.reader(csvfile, delimiter=',')
            for row in region_data:
                id = int(row[1])
                iso = str(row[2])
                new_location = Location(iso, (0.0, 0.0))
                locations = []
                locations.append(new_location)
                population_size = int(row[3])
                number_of_agents = max(int(population_size * self.scale_factor), 1)
                total_number_of_agents += number_of_agents
                agents = []
                for _ in range(number_of_agents):
                    new_agent = Agent(0)
                    new_agent.weekly_routine = ["Default" for _ in range(ticks_in_week)]
                    new_agent.activity_locations = {"Default": [new_location]}
                    new_agent.activity_location_weights = {"Default": [1.0]}
                    agents.append(new_agent)
                new_region = Region(id, iso, activities, agents, locations)
                world.regions.append(new_region)
        world.number_of_regions = len(world.regions)

        log.info("Created world with %d agents", total_number_of_agents)
