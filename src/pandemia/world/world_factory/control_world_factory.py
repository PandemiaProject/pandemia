"""This file creates a test world."""

import logging
import random
import numpy as np
from collections import defaultdict

from ..agent import Agent 
from ..location import Location
from ..region import Region
from .. import World
from . import WorldFactory

log = logging.getLogger('test_world_factory')

class ControlWorldFactory(WorldFactory):
    """Creates a world which can be used for testing all components of the model (individually or collectively)."""

    def __init__(self, config, clock, scale_factor):

        self.config = config
        self.clock = clock
        self.scale_factor = scale_factor

        self.contacts            = self.config['contacts']
        self.contact_hours       = self.config['contact_hours']
        self.travel_prob_per_day = self.config['travel_prob_per_day']

    def get_world(self) -> World:

        log.info("Creating world...")

        world = World(self.scale_factor)

        self._create_test_regions(world)

        world.travel_matrix =\
            self.get_travel_matrix(world, self.travel_prob_per_day)

        return world

    def _create_test_regions(self, world):
        """Creates test agents and test locations and assembles them into test regions"""

        random_seed                   = self.config['random_seed']
        region_names                  = self.config['region_names']
        number_of_regions             = len(region_names)

        # Agents
        max_age                       = self.config['max_age']
        number_of_agents_per_region   = self.config['number_of_agents_per_region']
        activities                    = self.config['activities']
        activity_location_types       = self.config['activity_location_types']
        number_of_locations_for_agent = self.config['number_of_locations_for_agent_by_activity']

        # Locations
        location_type_counts          = self.config['location_type_counts']

        # Rescale
        number_of_agents_per_region = max(int(number_of_agents_per_region * self.scale_factor), 1)
        for typ in location_type_counts:
            location_type_counts[typ] =\
                max(int(location_type_counts[typ] * self.scale_factor), 1)

        ticks_in_week = self.clock.ticks_in_week
        ticks_in_day = self.clock.ticks_in_day

        world.regions = []
        world.number_of_regions = number_of_regions

        random.seed(random_seed)

        total_number_of_agents = 0

        log.info("Creating %i test regions...", number_of_regions)
        for id in range(number_of_regions):

            # Create locations

            locations = []

            locations_by_type = defaultdict(list)

            for typ in location_type_counts:
                for _ in range(location_type_counts[typ]):
                    coord = (0.0, 0.0)
                    new_location = Location(typ, coord)
                    locations_by_type[typ].append(new_location)

            for typ in locations_by_type:
                locations += locations_by_type[typ]

            # Create agents

            agents = []

            valid_locations = defaultdict(list)

            for activity in activities:
                for typ in activity_location_types[activity]:
                    valid_locations[activity] += locations_by_type[typ]

            total_number_of_agents += number_of_agents_per_region

            for _ in range(number_of_agents_per_region):
                age = random.randrange(max_age + 1)
                new_agent = Agent(age)
                new_agent.weekly_routine =\
                    self._random_routine(activities, ticks_in_day, ticks_in_week)
                for activity in activities:
                    num = min(number_of_locations_for_agent[activity],
                              len(valid_locations[activity]))
                    random_sample = random.sample(valid_locations[activity], num)
                    new_agent.activity_locations[activity] = random_sample
                    new_agent.activity_location_weights[activity] = [1/num for _ in range(num)]
                agents.append(new_agent)

            name = region_names[id]

            world.regions.append(Region(id, name, activities, agents, locations))

        log.info("Created world with %d agents", total_number_of_agents)

    def _random_routine(self, activities, ticks_in_day, ticks_in_week):
        """Creates random weekly routine"""

        tick = 0
        test_routine = [random.choices(activities, k=1)[0]]
        while tick < ticks_in_week - 1:
            activity = random.choices(activities, k=1)[0]
            duration = min(random.randrange(ticks_in_week - tick), ticks_in_day)
            test_routine += [activity for _ in range(duration)]
            tick += duration

        return test_routine

    def get_travel_matrix(self, world, travel_prob_per_day):
        """Constructs matrix of mixing between regions."""

        regions = world.regions

        num_of_regions = len(regions)

        # This will be the matrix returned by the function
        baseline_agents_travelling_matrix =\
            np.zeros((num_of_regions, num_of_regions), dtype=np.int64)

        # Calculate local travel, rescaled according to step size
        share_matrix = np.zeros((num_of_regions, num_of_regions), dtype=np.float64)
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
