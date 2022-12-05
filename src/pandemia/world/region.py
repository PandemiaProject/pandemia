"""Represents a region consisting of activities, agents and locations"""

import logging
import numpy as np

from .agent import Agent
from .location import Location

log = logging.getLogger("region")

class VectorRegion:
    """Contains vector representations of agents, locations and other objects. Additional
    attributes are initialized by the simulation components, for example the health model,
    inside their vectorize_component function."""

    def __init__(self, id: int,
                       name: str,
                       ticks_in_week: int,
                       number_of_activities: int,
                       number_of_agents: int,
                       number_of_locations: int,
                       max_num_activity_locations: int):

        self.id = id
        self.name = name
        self.other_name = None
        self.super_region = None
        self.random_state = None
        self.prng = None
        self.number_of_agents = number_of_agents
        self.age = np.zeros((number_of_agents), dtype=np.int64)
        self.number_of_locations = number_of_locations

        assert number_of_activities <= 255
        self.number_of_activities = number_of_activities
        self.weekly_routines = np.zeros((number_of_agents, ticks_in_week), dtype=np.uint8)

        self.num_activity_locations =\
            np.zeros((number_of_agents, number_of_activities), dtype=np.int64)
        self.activity_locations = np.zeros((number_of_agents, number_of_activities,
                                            max_num_activity_locations), dtype=np.int64)
        self.activity_location_weights = np.zeros((number_of_agents, number_of_activities,
                                                   max_num_activity_locations), dtype=np.float64)
        self.max_num_activity_locations = max_num_activity_locations
        self.activity_strings = [None for _ in range(number_of_activities)]
        self.location_typ_strings = [None for _ in range(number_of_locations)]
        self.location_x_coords = np.zeros((number_of_locations), dtype=np.float64)
        self.location_y_coords = np.zeros((number_of_locations), dtype=np.float64)
        self.coordinates = None
        self.region_coordinates = None

        # These members are not consistently applied between Void and Default components
        self.current_disease = np.zeros((number_of_agents), dtype=float)
        self.current_strain = np.full((number_of_agents), -1, dtype=int)
        self.requested_location_update = np.zeros((number_of_agents), dtype=int)
        self.requesting_location_update = np.zeros((number_of_agents), dtype=int)
        self.home_location = np.zeros((number_of_agents), dtype=int)
        self.current_quarantine = np.zeros((number_of_agents), dtype=int)
        self.current_location = np.zeros((number_of_agents), dtype=int)
        self.current_facemask = np.zeros((number_of_agents), dtype=int)
        self.requested_facemask_update = np.zeros((number_of_agents), dtype=int)
        self.requesting_facemask_update = np.zeros((number_of_agents), dtype=int)
        self.wears_facemask = np.zeros((number_of_agents, number_of_activities), dtype=int)

class Region:
    """Represents a region"""

    def __init__(self, id: int, name: str, activities: list[str],
                 agents: list[Agent], locations: list[Location]):

        self.id = id
        self.name = name
        self.other_name = None
        self.super_region = None
        self.activities = activities
        self.agents = agents
        self.locations = locations
        self.coordinates = None
        self.region_coordinates = None

    def vectorize_region(self):
        """Converts object of type Region to object of type VectorRegion"""

        number_of_activities = len(self.activities)
        activity_ids = {}
        for a in range(len(self.activities)):
            activity_ids[self.activities[a]] = a

        number_of_agents = len(self.agents)
        agent_ids = {}
        for n in range(len(self.agents)):
            agent_ids[self.agents[n]] = n

        number_of_locations = len(self.locations)
        location_ids = {}
        for m in range(number_of_locations):
            location_ids[self.locations[m]] = m

        ticks_in_week = len(self.agents[0].weekly_routine)
        for agent in self.agents:
            assert ticks_in_week == len(agent.weekly_routine)

        max_num_activity_locations = 0
        for agent in self.agents:
            max_for_agent = max([len(agent.activity_locations[a]) for a in self.activities])
            if max_for_agent > max_num_activity_locations:
                max_num_activity_locations = max_for_agent

        vector_region = VectorRegion(self.id, self.name, ticks_in_week, number_of_activities,
                                     number_of_agents, number_of_locations,
                                     max_num_activity_locations)

        vector_region.other_name = self.other_name
        vector_region.super_region = self.super_region

        for agent in self.agents:
            agent_id = agent_ids[agent]
            vector_region.age[agent_id] = agent.age
            weekly_routine_ids = [activity_ids[activity] for activity in agent.weekly_routine]
            vector_region.weekly_routines[agent_id] = weekly_routine_ids
            for activity in self.activities:
                activity_id = activity_ids[activity]
                loc_ids = [location_ids[location] for location
                            in agent.activity_locations[activity]]
                vector_region.num_activity_locations[agent_id][activity_id] = len(loc_ids)
                weights = agent.activity_location_weights[activity]
                for j in range(len(loc_ids)):
                    vector_region.activity_locations[agent_id][activity_id][j] = loc_ids[j]
                    vector_region.activity_location_weights[agent_id][activity_id][j] = weights[j]
                for j in range(len(loc_ids), max_num_activity_locations):
                    vector_region.activity_locations[agent_id][activity_id][j] = -1
                    vector_region.activity_location_weights[agent_id][activity_id][j] = 0.0

        # Determine location type strings and coordinates
        for m in range(number_of_locations):
            vector_region.location_typ_strings[m] = self.locations[m].typ
            vector_region.location_x_coords[m] = self.locations[m].coord[0]
            vector_region.location_y_coords[m] = self.locations[m].coord[1]

        # Determine activity strings
        for a in range(number_of_activities):
            vector_region.activity_strings[a] = self.activities[a]

        # Determine shape data
        vector_region.coordinates = self.coordinates
        vector_region.region_coordinates = self.region_coordinates

        return vector_region
