"""Represents a region"""

import logging
import numpy as np

from .agent import Agent
from .location import Location

log = logging.getLogger("region")

class Region:

    def __init__(self, id: int, name: str, activities: list[str],
                 agents: list[Agent], locations: list[Location]):
        """Represents a region, for example a country or an administrative division, consisting of
        agents, locations and activities.

        Attributes:
            id (int): An integer identifier for this region.
            name (str): The name of the region. For example, if the region represents a country,
                this could be the country code in ISO 3166-1 alpha-2 format.
            other_name (Union[NoneType, str]): Another name of the region. For example, if the
                region represents a country, this could be the country code in ISO 3166-1 alpha-3
                format.
            super_region (Union[NoneType, str]): The group of regions to which this region belongs.
                For example, if the regions represents a country, its super region might be the
                continent to which it belongs.
            activities (list[str]): A list of all activities performed by agents in this region.
            agents (list[Agent]): A list of all agents in this region.
            locations (list[Location]): A list of all locations in this region.
            coordinates (Union[NoneType, list]): Used to render the region as polygons. If provided,
                the list should be of the format:
              
                    [[points_0], [points_1], ..., [points_N]]

                with each points_n a list of 2-tuples of floats. A list points_n represents the x, y
                coordinates of each point along the border of a connected piece of the region.
            region_coordinates (Union[NoneType, list]): Used to render the region as grid squares.
                If provided, the list should be a list
                of 2-tuples of floats, representing the x, y coordinates of each grid square.
        """

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
        """Converts an object of type Region to an object of type VectorRegion. The latter objects
        are vectorized versions of the former, consisting mainly of numpy arrays.
        
        Returns:
          out (VectorRegion): A vector representation of the region.
        """

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

        # Determine coordinates data
        vector_region.coordinates = self.coordinates
        vector_region.region_coordinates = self.region_coordinates

        return vector_region

class VectorRegion:

    def __init__(self, id: int,
                       name: str,
                       ticks_in_week: int,
                       number_of_activities: int,
                       number_of_agents: int,
                       number_of_locations: int,
                       max_num_activity_locations: int):

        """Represents a region, for example a country or an administrative division, in vector
        format. Agents, locations and activities are now represented as integers, with various
        arrays used to store their attributes. Further attributes, in addition to the ones below,
        may be added to the vector world by simulation components, where necessary.

        Attributes:
            id (int): An integer identifier for this region.
            name (str): The name of the region. For example, if the region represents a country,
                this could be the country code in ISO 3166-1 alpha-2 format.
            other_name (Union[NoneType, str]): Another name of the region. For example, if the
                region represents a country, this could be the country code in ISO 3166-1 alpha-3
                format.
            super_region (Union[NoneType, str]): The group of regions to which this region belongs.
                For example, if the regions represents a country, its super region might be the
                continent to which it belongs.
            random_state (Union[NoneType, numpy.ndarray]): A pair of 64-bit integers. This pair is
                the random seed for the prng used inside the C libraries. This is an numpy array of
                length 2, consisting of two integers of type numpy.uint64. Each region gets it own
                random state, to preserve determinism when parallelizing.
            prng (Random): An instance of the Random class for this region. Each region gets it own
                instance of the Random class, to preserve determinism when parallelizing.
            number_of_agents (int): The number of agents in this region.
            age (numpy.ndarray): An array of length number_of_agents recording the age of each
                agent.
            number_of_locations (int): The number of locations in this region.
            number_of_activities (int): The number of activities performed by agents in this region.
            weekly_routines (numpy.ndarray): For each agent, a sequence of integers of length
                ticks_in_week, indicating which activities (represented as integers) are performed
                each tick.
            num_activity_locations (numpy.ndarray): For each agent, for each activity, the number of
                locations at which that agent can perform that activitiy.
            activity_locations (numpy.ndarray): For each agent, for each activity, the locations
                (represented as integers) at which the agent can perform that activity. The number
                of locations at which agents perform activities may be variable, as recorded by
                num_activity_locations, but should be bounded above by max_num_activity_locations.
                This array is not a jagged array, the entries between num_activity_locations and
                max_num_activity_locations being padding.
            activity_location_weights (numpy.ndarray): For each agent, for each activity, the
                weights for each location at which the agent can perform that activity. The weights
                need not sum to one, with these weights being only relavent upto
                num_activity_locations. This array is not a jagged array, the entries between
                num_activity_locations and max_num_activity_locations being padding.
            max_num_activity_locations (int): The maximum number of locations at which it is
                possible to perform an activity. This maximum is taken accross all agents and all
                activities in this region.
            activity_strings (list[str]): The list of activities for this region. If an activity is
                represented by the integer a, then the corresponding string will be
                activity_strings[a], for example "Home".
            location_typ_strings (list[str]): If a location is represented by the integer l, then
                location_typ_strings[l] gives the type of location l, for example "House".
            location_x_coords (numpy.ndarray): The x coordinates of all locations in this region.
            location_y_coords (numpy.ndarray): The y coordinates of all locations in this region.
            coordinates (Union[NoneType, list]): Used to render the region as polygons. If provided,
                the list should be of the format
              
                    [[points_0], [points_1], ..., [points_N]]

                with each points_n a list of 2-tuples of floats. A list points_n represents the x, y
                coordinates of each point along the border of a connected piece of the region.
            region_coordinates (Union[NoneType, list]): Used to render the region as grid squares.
                If provided, the list should be a list of 2-tuples of floats, representing the x, y
                coordinates of each grid square.
        """

        self.id = id
        self.name = name
        self.other_name = None
        self.super_region = None
        self.random_state = None
        self.prng = None
        self.number_of_agents = number_of_agents
        self.age = np.zeros((number_of_agents), dtype=int)
        self.number_of_locations = number_of_locations

        assert number_of_activities <= 255
        self.number_of_activities = number_of_activities
        self.weekly_routines = np.zeros((number_of_agents, ticks_in_week), dtype=np.uint8)

        self.num_activity_locations =\
            np.zeros((number_of_agents, number_of_activities), dtype=int)
        self.activity_locations = np.zeros((number_of_agents, number_of_activities,
                                            max_num_activity_locations), dtype=int)
        self.activity_location_weights = np.zeros((number_of_agents, number_of_activities,
                                                   max_num_activity_locations), dtype=float)
        self.max_num_activity_locations = max_num_activity_locations
        self.activity_strings = [None for _ in range(number_of_activities)]
        self.location_typ_strings = [None for _ in range(number_of_locations)]
        self.location_x_coords = np.zeros((number_of_locations), dtype=float)
        self.location_y_coords = np.zeros((number_of_locations), dtype=float)
        self.coordinates = None
        self.region_coordinates = None
