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

class ColombiaWorldFactory(WorldFactory):
    """Reads a DensityMap and generates a world based on the densities indicated therein."""

    def __init__(self, config, clock, scale_factor):
        """Create agents and locations according to the population density map given"""

        self.config = config
        self.clock = clock
        self.scale_factor = scale_factor

        self.regions_shape_path        = self.config['regions_shape_path']
        self.regions_data_path         = self.config['regions_data_path']
        self.local_travel_prob_per_day = self.config['local_travel_prob_per_day']

    def get_world(self) -> World:

        log.info("Creating world...")

        world = World(self.scale_factor)

        self._create_test_regions(world)

        world.travel_matrix =\
            self.get_travel_matrix(world, self.local_travel_prob_per_day)

        return world

    def _create_test_regions(self, world):
        """Creates world for example one"""

        ticks_in_week = self.clock.ticks_in_week

        # Initialize regions
        self.initialize_regions(self.regions_data_path, ticks_in_week, world)

        # Add shape data to regions
        self.add_shape_data(self.regions_shape_path, world)

    def initialize_regions(self, regions_data_path, ticks_in_week, world):
        """Initialize regions"""

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

    def add_shape_data(self, regions_shape_path, world):
        """Add polgonal shapes to regions for rendering"""

        sf = shp.Reader(regions_shape_path)
        shape_recs = sf.shapeRecords()

        # Determine for which regions coordinate data can be found
        regions_to_coordinates = {}
        for region in world.regions:
            regions_to_coordinates[region] = None
            for id in range(len(shape_recs)):
                iso = shape_recs[id].record[1]
                if region.name == iso:
                    # Extract coordinates according to shape type
                    geoj = sf.shape(id).__geo_interface__
                    if geoj["type"] == "Polygon":
                        coordinates = [geoj["coordinates"]]
                    elif geoj["type"] == "MultiPolygon":
                        coordinates = geoj["coordinates"]
                    else:
                        coordinates = None
                    region.coordinates = coordinates

    def get_travel_matrix(self, world, travel_prob_per_day):
        """Constructs matrix of mixing between regions."""

        regions = world.regions

        num_of_regions = len(regions)

        # This will be the matrix returned by the function
        baseline_agents_travelling_matrix =\
            np.zeros((num_of_regions, num_of_regions), dtype=int)

        # Calculate local travel, rescaled according to step size
        share_matrix = np.zeros((num_of_regions, num_of_regions), dtype=float)
        for region in regions:
            for other_region in regions:
                if region.id != other_region.id:
                    share_matrix[region.id][other_region.id] = len(other_region.agents)
        for i in range(num_of_regions):
            row_sum = sum(share_matrix[i])
            for j in range(num_of_regions):
                if share_matrix[i][j] > 0:
                    share_matrix[i][j] = share_matrix[i][j] / row_sum
        ids_to_population_sizes = {region.id: len(region.agents) for region in regions}

        for i in range(num_of_regions):
            travel_per_day = ids_to_population_sizes[i] * travel_prob_per_day
            for j in range(num_of_regions):
                baseline_agents_travelling_matrix[i][j] = share_matrix[i][j] * travel_per_day

        return baseline_agents_travelling_matrix
