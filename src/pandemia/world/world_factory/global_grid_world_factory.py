"""This file creates a world."""

import logging
import numpy as np
import os
import csv

from collections import defaultdict
from tqdm import tqdm
from scipy.spatial import KDTree
import shapefile as shp
from scipy.spatial import distance
import multiprocessing
from ..agent import Agent 
from ..location import Location
from ..region import Region
from .. import World
from . import WorldFactory
from ...random_tools import Random

log = logging.getLogger('global_grid_world_factory')

class GlobalGridWorldFactory(WorldFactory):
    """In this World model, for each region there are two activities, namely Home_Activity and
    Community_Activity. Each agent in assigned a weekly routine, with days divided into three
    periods, each period lasting eight hours. Each agent is assigned a home, at which they perform
    the activity Home_Activity, with the homes spatially distributed within each region using
    population grid data for that region. Each agent is also assigned a set of grid squares, which
    includes the square in which their assigned home is located. When performing the activity
    Community_Activity, an agent will choose from these grid squares. The grid squares are assigned
    to each agent and weighted in such a way that agents are more likely to visit squares close to
    their home that squares far away. The weights are chosen according to the so-called gravity
    model of human mobility. There are therefore two types of location in this model, namely
    House and Square. Locations of type Square are assigned a location specific transmission
    multiplier. This is a free parameter that modulates the intensity of community transmission. In
    this model, travel between regions is based on air travel data.
    """

    def __init__(self, config, clock, scale_factor):

        self.config = config
        self.clock = clock
        self.scale_factor = scale_factor

        # Random seed
        self.random_seed                 = self.config['random_seed']

        # Region to include in the simulation
        self.regions_to_simulate         = self.config['regions_to_simulate']

        # Data files
        self.regions_data_file           = self.config['regions_data_file']
        self.regions_shape_data_file     = self.config['regions_shape_data_file']
        self.population_dist_data_folder = self.config['population_dist_data_folder']
        self.population_id_data_folder   = self.config['population_id_data_folder']
        self.airport_data_file           = self.config['airport_data_file']
        self.air_travel_data_file        = self.config['air_travel_data_file']

        # Agents, locations and activities
        self.house_location_type         = self.config['house_location_type']
        self.square_location_type        = self.config['square_location_type']
        self.home_activity               = self.config['home_activity']
        self.community_activity          = self.config['community_activity']
        self.min_age_community_activity  = self.config['min_age_community_activity']
        self.max_age_community_activity  = self.config['max_age_community_activity']
        self.local_sample_size           = self.config['local_sample_size']
        self.nonlocal_sample_size        = self.config['nonlocal_sample_size']
        self.subsample_size              = self.config['num_locs_for_community_activity']
        self.gamma                       = self.config['gravity_model_exponent']

        # Travel matrix
        self.local_travel_prob_per_day   = self.config['local_travel_prob_per_day']
        self.distance_threshold          = self.config['distance_threshold']

    def get_world(self) -> World:

        log.info("Creating world...")

        world = World(self.scale_factor)

        # Create regions using population data
        self.region_data = {}
        id = 0
        with open(self.regions_data_file, newline='') as csvfile:
            next(csvfile)
            region_data = csv.reader(csvfile, delimiter=',')
            for row in region_data:
                iso2 = str(row[1])
                if (len(self.regions_to_simulate) == 0) or (iso2 in self.regions_to_simulate):
                    iso3 = str(row[2])
                    super_region = str(row[3])
                    household_size = round(float(row[4]))
                    population_size = int(row[5])
                    age_distribution = [float(row[6 + r]) for r in range(101)]
                    number_of_agents = max(int(population_size * self.scale_factor), 1)
                    new_region = self._create_region(id, iso2, iso3, super_region, number_of_agents,
                                                        age_distribution, household_size)
                    world.regions.append(new_region)
                    id += 1

        world.number_of_regions = len(world.regions)

        number_of_agents = sum([len(region.agents) for region in world.regions])

        log.info("Created world with %d agents", number_of_agents)

        self.add_shape_data(self.regions_shape_data_file, world)

        world.travel_matrix =\
            self.get_travel_matrix(world, self.airport_data_file, self.air_travel_data_file,
                                   self.local_travel_prob_per_day,
                                   self.distance_threshold)

        return world

    def _create_region(self, id, iso2, iso3, super_region, number_of_agents,
                       age_distribution, household_size):
        """Creates test agents and test locations and assembles them into test regions"""

        assert self.clock.ticks_in_day == 3

        log.info("Creating region: " + iso2 + "...")

        prng = Random(self.random_seed)

        log.debug("Loading map...")

        # Extract coordinates of top left corner
        for filename in os.listdir(self.population_dist_data_folder):
            if iso3 == filename[0:3].upper():
                filepath_dist = os.path.join(self.population_dist_data_folder, filename)
                if os.path.isfile(filepath_dist):
                    with open(filepath_dist, newline='') as csvfile:
                        region_data = csv.reader(csvfile, delimiter=' ')
                        ncols = int(next(region_data)[-1])
                        nrows = int(next(region_data)[-1])
                        xll = int(next(region_data)[-1])
                        yll = int(next(region_data)[-1])
                        square_size = float(next(region_data)[-1])
                        xtl = xll
                        ytl = yll + (square_size * nrows)
        for filename in os.listdir(self.population_id_data_folder):
            if iso3 == filename[0:3].upper():
                filepath_id = os.path.join(self.population_id_data_folder, filename)

        # Load population distribution map
        map_dist = np.genfromtxt(filepath_dist, delimiter=' ', skip_header=6, dtype=float)
        map_dist /= np.sum(map_dist)
        map_id = np.genfromtxt(filepath_id, delimiter=' ', skip_header=6, dtype=float)

        # Determine which squares belong to each region, for plotting purposes
        region_x_coords = []
        region_y_coords = []
        for y in range(map_id.shape[0]):
            for x in range(map_id.shape[1]):
                if map_id[y][x] > 0:
                    region_x_coords.append(xtl + (square_size * x))
                    region_y_coords.append(ytl - (square_size * y))
        region_coordinates = zip(region_x_coords, region_y_coords)

        # The ideal (float) and rounded (int) distributions of population
        map_floats = number_of_agents * map_dist
        map_ints = map_floats.astype(int)

        # We must account for the rounding error in order to generate the correct number of agents.
        # The first step in accounting for this error involves adding on multiples of the map_dist.
        log.debug("Distributing population (approximate)...")
        old_sum = -1
        converged = False
        while ((np.sum(map_ints) < number_of_agents) and not converged):
            err = np.sum(map_floats - map_ints)
            addition = err * map_dist
            map_ints = (map_ints + addition).astype(int)
            new_sum = np.sum(map_ints)
            if new_sum == old_sum:
                converged = True
            old_sum = new_sum

        # The second step involves adding on single agents in such a way that minimizes the error.
        # Doing this step alone is far too slow, hence the previous step.
        log.debug("Distributing population (exact)...")
        while (np.sum(map_ints) < number_of_agents):
            err = np.square(map_ints - map_floats)
            map_ints[np.unravel_index(np.argmax(err), err.shape)] += 1
        self.map = map_ints

        # Check that the correct number of agents will be created via this process
        assert np.sum(self.map) == number_of_agents

        # Create activities

        activities = [self.home_activity, self.community_activity]

        # Create grid squares and houses and populate the houses with agents, assigning them a
        # location for the home activity

        log.debug("Creating locations and agents...")

        locations = []
        agents = []
        squares = []
        self.squares_to_agents = defaultdict(list)

        for y in range(self.map.shape[0]):
            for x in range(self.map.shape[1]):
                if self.map[y][x] > 0:

                    # Create grid square
                    new_square = Location(self.square_location_type, (xtl + (square_size * x),
                                                                      ytl - (square_size * y)))
                    locations.append(new_square)
                    squares.append(new_square)

                    # Create agents to live in this square
                    unhoused = []
                    for _ in range(self.map[y][x]):
                        age = prng.random_choices(list(range(101)), age_distribution, 1)[0]
                        new_agent = Agent(age)
                        self.squares_to_agents[new_square].append(new_agent)
                        new_agent.weekly_routine = self._assign_routine(age)
                        agents.append(new_agent)
                        unhoused.append(new_agent)

                    # Partition these agents into households and create houses
                    num_houses_needed = (self.map[y][x] // household_size) + 1
                    prng.random_shuffle(unhoused)
                    households = [unhoused[n::num_houses_needed] for n in range(num_houses_needed)]
                    for household in households:
                        new_house = Location(self.house_location_type, (xtl + (square_size * x),
                                                                        ytl - (square_size * y)))
                        locations.append(new_house)
                        for agent in household:
                            agent.activity_locations[self.home_activity] = [new_house]
                            agent.activity_location_weights[self.home_activity] = [1.0]

        # Assign agents locations for the community activity based on weighted proximity

        log.debug("Assigning community locations...")

        kdtree = KDTree([square.coord for square in squares])
        local_sample_size = min(self.local_sample_size, len(squares))
        nonlocal_sample_size = min(self.nonlocal_sample_size, len(squares))

        for square in tqdm(squares):

            if local_sample_size > 1:
                _, nearest_indices = kdtree.query(square.coord, local_sample_size)
                local_sample = [squares[i] for i in nearest_indices]
            else:
                local_sample = [square]

            nonlocal_sample = prng.random_sample(squares, nonlocal_sample_size)

            sample = local_sample + nonlocal_sample

            weights = []
            for other_square in sample:
                weights.append(self._travel_weight(square, other_square))

            activity = self.community_activity

            for agent in self.squares_to_agents[square]:

                # All agents can perform the community activity in their home square
                agent.activity_locations[activity] = [square]
                agent.activity_location_weights[activity] =\
                    [self._travel_weight(square, square)]

                # In addition to their home square, other squares are added by weighted sampling
                if self.subsample_size > 1:
                    subsample_indexes =\
                        prng.random_sample(range(len(sample)),
                                           min(self.subsample_size - 1, len(sample)))
                    subsample_locations =\
                        [sample[index] for index in subsample_indexes]
                    subsample_weights =\
                        [weights[index] for index in subsample_indexes]
                    agent.activity_locations[activity].extend(subsample_locations)
                    agent.activity_location_weights[activity].extend(subsample_weights)

        new_region = Region(id, iso2, activities, agents, locations)

        new_region.other_name = iso3
        new_region.super_region = super_region

        new_region.region_coordinates = region_coordinates

        return new_region

    def _travel_weight(self, square, other_square):
        """For a given distance, returns a weight"""

        dist = square.distance_euclidean(other_square)
        population_other_square = len(self.squares_to_agents[other_square])

        weight = population_other_square / (1 + dist**self.gamma)

        return weight

    def _assign_routine(self, age):
        """Creates random weekly routine"""

        weekend_day = [self.home_activity, self.home_activity, self.home_activity]
        week_day = [self.home_activity, self.community_activity, self.home_activity]

        if age >= self.min_age_community_activity and age <= self.max_age_community_activity:

            weekly_routine = weekend_day + (week_day * 5) + weekend_day

        else:

            weekly_routine = weekend_day * 7

        return weekly_routine

    def add_shape_data(self, regions_shape_data_file, world):
        """Add polgonal shapes to regions for rendering"""

        sf = shp.Reader(regions_shape_data_file)
        shape_recs = sf.shapeRecords()

        # Determine for which regions coordinate data can be found
        for region in world.regions:
            for id in range(len(shape_recs)):
                iso = shape_recs[id].record[0]
                if iso == 'EL':
                    iso = 'GR'
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

    def get_travel_matrix(self, world, airport_data_file, air_travel_data_file,
                          local_travel_prob_per_day, distance_threshold):
        """Constructs matrix of mixing between regions, constructed via a combination of
        air travel and local travel. Since the air travel data used for these simulaitons also
        records the month of travel, we additionally calculate an air travel matrix for each month,
        for later use. Note that the matrix gets rescaled, using a scale_factor, not here but in
        the regional mixing model."""

        regions = world.regions

        num_of_regions = len(regions)
        months_in_year = 12

        # This will be the matrix returned by the function
        baseline_agents_travelling_matrix =\
            np.zeros((num_of_regions, num_of_regions), dtype=float)

        # It will be constructed by obtaining air travel with local travel, the latter meaning
        # travel to neighbouring regions
        daily_air_travel_matrix = np.zeros((num_of_regions, num_of_regions), dtype=float)
        daily_air_travel_matrix_by_month =\
            np.zeros((months_in_year, num_of_regions, num_of_regions), dtype=float)
        daily_local_travel_matrix = np.zeros((num_of_regions, num_of_regions), dtype=float)

        # Map airports to region isos
        airports_to_region_iso = {}
        with open(airport_data_file, newline='') as csvfile:
            airport_data = csv.reader(csvfile, delimiter=',')
            next(airport_data, None)
            for row in airport_data:
                airport = str(row[0])
                region_iso = str(row[3])
                airports_to_region_iso[airport] = region_iso

        # Get air travel matrix, recording the number of travellers between regions either per month
        # or year
        region_isos = [region.name for region in regions]
        region_isos_to_ids = {region.name: region.id for region in regions}
        air_travel_matrix_by_month =\
            np.zeros((months_in_year, num_of_regions, num_of_regions), dtype=int)
        air_travel_matrix = np.zeros((num_of_regions, num_of_regions), dtype=int)
        with open(air_travel_data_file, newline='') as csvfile:
            travel_data = csv.reader(csvfile, delimiter=',')
            next(travel_data, None)
            for row in travel_data:
                origin_airport = str(row[0])
                destination_airport = str(row[1])
                origin_region_iso = airports_to_region_iso[origin_airport]
                destination_region_iso = airports_to_region_iso[destination_airport]
                if (origin_region_iso in region_isos) and (destination_region_iso in region_isos):
                    origin_region_id = region_isos_to_ids[origin_region_iso]
                    destination_region_id = region_isos_to_ids[destination_region_iso]
                    # Discard internal flights
                    if origin_region_id != destination_region_id:
                        m = int(row[2]) - 1
                        prediction = int(row[3])
                        air_travel_matrix_by_month[m][origin_region_id][destination_region_id] +=\
                            prediction
                        air_travel_matrix[origin_region_id][destination_region_id] += prediction

        # Calculate air travel, rescaled from months or years to days
        days_in_year_2010 = 365
        days_in_month_2010 = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        for month in range(months_in_year):
            for i in range(num_of_regions):
                for j in range(num_of_regions):
                    air_travel_per_day =\
                        air_travel_matrix_by_month[month][i][j] / days_in_month_2010[month]
                    daily_air_travel_matrix_by_month[month][i][j] = air_travel_per_day
        for month in range(months_in_year):
            for i in range(num_of_regions):
                daily_air_travel_matrix_by_month[month][i][i] = 0
        for i in range(num_of_regions):
            for j in range(num_of_regions):
                air_travel_per_day = air_travel_matrix[i][j] / days_in_year_2010
                daily_air_travel_matrix[i][j] = air_travel_per_day
        np.fill_diagonal(daily_air_travel_matrix, 0)

        # Get adjacency matrix, recording which regions border which others
        adjacency_matrix = np.zeros((num_of_regions, num_of_regions), dtype=int)
        for region in regions:
            for other_region in regions:
                if (region.coordinates is not None) and (other_region.coordinates is not None):
                    list1 = []
                    for part in region.coordinates:
                        list1 += part[0]
                    list2 = []
                    for part in other_region.coordinates:
                        list2 += part[0]
                    nplist1 = np.array(list1)
                    nplist2 = np.array(list2)
                    dist = min(distance.cdist(nplist1,nplist2).min(axis=1))
                    if dist < distance_threshold:
                        adjacency_matrix[region.id][other_region.id] = 1
        np.fill_diagonal(adjacency_matrix, 0)

        # Calculate local travel, rescaled according to step size
        share_matrix = np.zeros((num_of_regions, num_of_regions), dtype=float)
        ids_to_population_sizes =\
            {region.id: len(region.agents) * (1 / self.scale_factor)\
                        for region in regions}
        for region in regions:
            for other_region in regions:
                if adjacency_matrix[region.id][other_region.id] == 1:
                    share_matrix[region.id][other_region.id] =\
                        ids_to_population_sizes[other_region.id]
        for i in range(num_of_regions):
            row_sum = sum(share_matrix[i])
            for j in range(num_of_regions):
                if share_matrix[i][j] > 0:
                    share_matrix[i][j] = share_matrix[i][j] / row_sum
        for i in range(num_of_regions):
            local_travel_per_day = ids_to_population_sizes[i] * local_travel_prob_per_day
            for j in range(num_of_regions):
                daily_local_travel_matrix[i][j] = share_matrix[i][j] * local_travel_per_day

        # Combine air travel and local travel probability matrices to get final travel
        for i in range(num_of_regions):
            for j in range(num_of_regions):
                baseline_agents_travelling_matrix[i][j] =\
                    daily_local_travel_matrix[i][j] + daily_air_travel_matrix[i][j]

        # Check there are not too many travellers
        for i in range(num_of_regions):
            assert sum(baseline_agents_travelling_matrix[i]) < ids_to_population_sizes[i]

        return baseline_agents_travelling_matrix
