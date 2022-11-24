"""This file creates a world."""

import logging
import csv
import numpy as np
import pickle

from ..agent import Agent 
from ..location import Location
from ..region import Region
from .. import World
from . import WorldFactory

log = logging.getLogger('world_factory')

class ABMluxWorldFactory(WorldFactory):
    """Takes a World in ABMlux format and converts it into Pandemia format. Worlds built using
    ABMlux can therefore serve as the basis for a Pandemia simulation."""

    def __init__(self, config, clock, scale_factor):

        self.config = config
        self.clock = clock
        self.scale_factor = scale_factor

        self.world_fp = self.config['world_fp']
        self.region_name = self.config['region_name']

    def get_world(self) -> World:

        log.info("Creating world...")

        world = World(self.scale_factor)

        self._create_region(world)

        number_of_regions = len(world.regions)

        assert number_of_regions == 1

        world.travel_matrix =\
            np.zeros((number_of_regions, number_of_regions), dtype=int)

        return world

    def _create_region(self, world):

        with open(self.world_fp, 'rb') as fin:
            ms_abmlux_sim_factory = pickle.load(fin)

        assert ms_abmlux_sim_factory.world.scale_factor == self.scale_factor
        assert ms_abmlux_sim_factory.clock.tick_length_s == self.clock.tick_length_s

        total_number_of_agents = len(ms_abmlux_sim_factory.world.agents)

        world.regions = []

        activities = list(ms_abmlux_sim_factory.activity_manager.map_config.keys())

        locations = []
        ms_abmlux_locs_to_pandemia_locs = {}
        for loc in ms_abmlux_sim_factory.world.locations:
            location_type = loc.typ
            coordinates = loc.coord
            new_location = Location(location_type, (coordinates[0], coordinates[1]))
            ms_abmlux_locs_to_pandemia_locs[loc] = new_location
            locations.append(new_location)

        agents = []
        for agent in ms_abmlux_sim_factory.world.agents:
            age = agent.age
            weekly_routine_int =\
                ms_abmlux_sim_factory.activity_model.weeks_by_agent[agent].weekly_routine
            weekly_routine_str =\
                [ms_abmlux_sim_factory.activity_manager.as_str(i) for i in weekly_routine_int]

            new_agent = Agent(age)
            new_agent.weekly_routine = weekly_routine_str

            for act_int in agent.activity_locations:
                abmlux_locs = agent.activity_locations[act_int]
                pandemia_locs = [ms_abmlux_locs_to_pandemia_locs[loc] for loc in abmlux_locs]
                act_str = ms_abmlux_sim_factory.activity_manager.as_str(act_int)
                new_agent.activity_locations[act_str] = pandemia_locs
                num_locs = len(pandemia_locs)
                new_agent.activity_location_weights[act_str] = [1/num_locs for _ in range(num_locs)]

            agents.append(new_agent)

        new_region = Region(0, self.region_name, activities, agents, locations)
        world.regions.append(new_region)
        world.number_of_regions = len(world.regions)

        log.info("Created world with %d agents", total_number_of_agents)
