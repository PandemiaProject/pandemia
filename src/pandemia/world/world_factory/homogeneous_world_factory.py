"""This file creates a world."""

import logging
import csv
import numpy as np
import shapefile as shp
from scipy.spatial import distance

from ..agent import Agent 
from ..location import Location
from ..region import Region
from .. import World
from . import WorldFactory

log = logging.getLogger('world_factory')

class HomogeneousWorldFactory(WorldFactory):
    """In this model, in each region there is only one activity, namely Default. There is only one
    location in each region. Agents mix within regions all within the same location. Using the
    sir_rescaling option in the default health model, agents in this world will mix homogeneously
    within each region. Travel between regions is determined using air travel data."""

    def __init__(self, config, clock, scale_factor):

        self.config = config
        self.clock = clock
        self.scale_factor = scale_factor

        self.regions_shape_path        = self.config['regions_shape_path']
        self.regions_data_path         = self.config['regions_data_path']
        self.airport_path              = self.config['airport_path']
        self.air_travel_path           = self.config['air_travel_path']
        self.local_travel_prob_per_day = self.config['local_travel_prob_per_day']
        self.distance_threshold        = self.config['distance_threshold']

    def get_world(self) -> World:

        log.info("Creating world...")

        world = World(self.scale_factor)

        self._create_test_regions(world)

        world.travel_matrix =\
            self.get_travel_matrix(world, self.airport_path, self.air_travel_path,
                                   self.local_travel_prob_per_day,
                                   self.distance_threshold)

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
        id = 0
        with open(regions_data_path, newline='') as csvfile:
            next(csvfile)
            region_data = csv.reader(csvfile, delimiter=',')
            for row in region_data:
                iso2 = str(row[1])
                iso3 = str(row[2])
                new_location = Location(iso2, (0.0, 0.0))
                locations = []
                locations.append(new_location)
                super_region = str(row[3])
                population_size = int(row[5])
                age_distribution = [float(row[6 + r]) for r in range(101)]
                number_of_agents = max(int(population_size * self.scale_factor), 1)
                total_number_of_agents += number_of_agents
                agents = []
                for age in range(101):
                    number_of_agents_of_this_age =\
                        max(int(age_distribution[age] * number_of_agents), 1)
                    for _ in range(number_of_agents_of_this_age):
                        new_agent = Agent(age)
                        new_agent.weekly_routine = ["Default" for _ in range(ticks_in_week)]
                        new_agent.activity_locations = {"Default": [new_location]}
                        new_agent.activity_location_weights = {"Default": [1.0]}
                        agents.append(new_agent)
                new_region = Region(id, iso2, activities, agents, locations)
                new_region.super_region = super_region
                new_region.other_name = iso3
                world.regions.append(new_region)
                id += 1
        world.number_of_regions = len(world.regions)

        log.info("Created world with %d agents", total_number_of_agents)

    def add_shape_data(self, regions_shape_path, world):
        """Add polgonal shapes to regions for rendering"""

        sf = shp.Reader(regions_shape_path)
        shape_recs = sf.shapeRecords()

        # Determine for which regions coordinate data can be found
        for region in world.regions:
            for id in range(len(shape_recs)):
                iso = shape_recs[id].record[0]
                if iso == 'UK':
                    iso = 'GB'
                if iso == 'RS':
                    iso = 'CS'
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
            np.zeros((num_of_regions, num_of_regions), dtype=np.float64)

        # It will be constructed by obtaining air travel with local travel, the latter meaning
        # travel to neighbouring regions
        daily_air_travel_matrix = np.zeros((num_of_regions, num_of_regions), dtype=np.float64)
        daily_air_travel_matrix_by_month =\
            np.zeros((months_in_year, num_of_regions, num_of_regions), dtype=np.float64)
        daily_local_travel_matrix = np.zeros((num_of_regions, num_of_regions), dtype=np.float64)

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
            np.zeros((months_in_year, num_of_regions, num_of_regions), dtype=np.int64)
        air_travel_matrix = np.zeros((num_of_regions, num_of_regions), dtype=np.int64)
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
        adjacency_matrix = np.zeros((num_of_regions, num_of_regions), dtype=np.int64)
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
        share_matrix = np.zeros((num_of_regions, num_of_regions), dtype=np.float64)
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
