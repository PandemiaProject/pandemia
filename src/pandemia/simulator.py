"""Simulates an epidemic"""

from collections import defaultdict
import logging
import numpy as np
import os
import csv
from datetime import datetime
from ctypes import c_int64, c_void_p, cdll
from joblib import Parallel, delayed
import platform
ext=".dll" if platform.system() == 'Windows' else ".so"
from .random_tools import Random

import uuid

from .version import VERSION

log = logging.getLogger('sim')

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class Simulator:
    """Class that simulates an outbreak."""

    def __init__(self,
                 config,
                 clock,
                 vector_world,
                 seasonal_effects_model,
                 health_model,
                 movement_model,
                 hospitalization_and_death_model,
                 testing_and_contact_tracing_model,
                 vaccination_model,
                 travel_model,
                 policy_maker_model,
                 telemetry_bus):

        # Static info
        self.pandemia_version = VERSION
        self.created_at       = datetime.now()
        self.run_id           = uuid.uuid4().hex
        self.lib = cdll.LoadLibrary(os.path.split(__file__)[0]+"/C/build/simulator_functions"+ext)

        log.info("Simulation created at %s with ID=%s", self.created_at, self.run_id)

        # Config
        self.config = config

        # Random seeds for each region
        self.random_seed = self.config['random_seed']

        # Telemetry bus used for reporting
        self.telemetry_bus = telemetry_bus

        # Components
        self.clock = clock
        self.vector_world = vector_world
        self.vector_regions = vector_world.vector_regions
        self.seasonal_effects_model = seasonal_effects_model
        self.health_model = health_model
        self.movement_model = movement_model
        self.hospitalization_and_death_model = hospitalization_and_death_model
        self.testing_and_contact_tracing_model = testing_and_contact_tracing_model
        self.vaccination_model = vaccination_model
        self.travel_model = travel_model
        self.policy_maker_model = policy_maker_model

        self.number_of_strains = self.health_model.number_of_strains
        self.number_of_vaccines = self.vaccination_model.number_of_vaccines
        self.number_of_regions = vector_world.number_of_regions

        self.collect_telemetry_data = self.lib.collect_telemetry_data
        self.count_dead             = self.lib.count_dead

        self.collect_telemetry_data.restype = None
        self.count_dead.restype             = c_int64

        # Parallel processing
        self.enable_parallel       = self.config['enable_parallel']
        self.num_jobs              = self.config['num_jobs']
        self.vector_region_batches = None

        # Output
        self.error = 0
        self.historical_data_filepath = self.config['historical_data_filepath']
        self.average_deaths_dict_historical = None
   
    def _set_telemetry_bus(self):
        """Assigns telemetry bus to each component"""

        self.seasonal_effects_model.set_telemetry_bus(self.telemetry_bus)
        self.health_model.set_telemetry_bus(self.telemetry_bus)
        self.movement_model.set_telemetry_bus(self.telemetry_bus)
        self.hospitalization_and_death_model.set_telemetry_bus(self.telemetry_bus)
        self.testing_and_contact_tracing_model.set_telemetry_bus(self.telemetry_bus)
        self.vaccination_model.set_telemetry_bus(self.telemetry_bus)
        self.travel_model.set_telemetry_bus(self.telemetry_bus)
        self.policy_maker_model.set_telemetry_bus(self.telemetry_bus)

    def _vectorize_components(self):
        """Initialize the numpy arrays associated to each component"""

        for vector_region in self.vector_regions:

            self.travel_model.vectorize_component(vector_region)
            self.policy_maker_model.vectorize_component(vector_region)
            self.seasonal_effects_model.vectorize_component(vector_region)
            self.health_model.vectorize_component(vector_region)
            self.movement_model.vectorize_component(vector_region)
            self.hospitalization_and_death_model.vectorize_component(vector_region)
            self.testing_and_contact_tracing_model.vectorize_component(vector_region)
            self.vaccination_model.vectorize_component(vector_region)

    def _parallelize_vector_regions(self):
        """Calculates number of jobs and partitions vector regions into approximately equally sized
        batches for parallel processing"""

        number_of_cpus = os.cpu_count()

        if self.num_jobs < 1:
            self.num_jobs = number_of_cpus
        elif self.num_jobs >= 1:
            self.num_jobs = min(self.num_jobs, number_of_cpus)
        else:
            self.num_jobs = 1

        max_agents   = max([v.number_of_agents for v in self.vector_regions])
        total_agents = sum([v.number_of_agents for v in self.vector_regions])

        ratio = (total_agents // max_agents) + 1

        self.num_jobs = min(self.num_jobs, self.number_of_regions, ratio)
        self.vector_region_batches = [[] for _ in range(self.num_jobs)]

        # Greedy partition algorithm
        for vector_region in sorted(self.vector_regions,
                                    key=lambda x: x.number_of_agents,
                                    reverse=True):
            self.vector_region_batches.sort(key=lambda x: sum([y.number_of_agents for y in x]))
            self.vector_region_batches[0].append(vector_region)

    def _seed_regions(self):
        """Sets random seeds for each region"""

        ULONG_MAX = 18446744073709551615

        # Create a prng for each vector region
        for vector_region in self.vector_regions:
            seed = self.random_seed
            # Each vector region gets an instance of Random()
            vector_region.prng = Random(seed)
            # Each vector region also gets a state for the prng used inside C functions
            np.random.seed(seed)
            vector_region.random_state =\
                np.random.randint(0, ULONG_MAX + 1, size=2, dtype=np.uint64)

    def _initial_conditions(self, offset):
        """Initialize the various submodels"""

        self.travel_model.initial_conditions(self)

        for vector_region in self.vector_regions:

            self.policy_maker_model.initial_conditions(vector_region)
            self.seasonal_effects_model.initial_conditions(vector_region)
            self.health_model.initial_conditions(vector_region)
            self.movement_model.initial_conditions(vector_region, offset)
            self.hospitalization_and_death_model.initial_conditions(vector_region)

            self._update(vector_region, 0)

        for vector_region in self.vector_regions:

            self.testing_and_contact_tracing_model.initial_conditions(vector_region)
            self.vaccination_model.initial_conditions(vector_region)

    def simulate_day(self, vector_regions, day, offset, ticks_in_day, ticks_in_week):
        """Simulate a day inside the given region"""

        for vector_region in vector_regions:

            self.seasonal_effects_model.dynamics(vector_region, day)
            self.policy_maker_model.dynamics(vector_region, day)

            for tick in range(ticks_in_day):

                t = (ticks_in_day * day) + tick

                self.health_model.dynamics(vector_region, t)
                self.movement_model.dynamics(vector_region, t, offset, ticks_in_week)
                self.hospitalization_and_death_model.dynamics(vector_region)

                self._update(vector_region, t)

            self.testing_and_contact_tracing_model.dynamics(vector_region, day)
            self.vaccination_model.dynamics(vector_region, day, ticks_in_day)

    def setup(self):
        """Setup a new simulation"""

        # Set the correct time
        self.clock.reset()
        offset = self.clock.ticks_through_week()

        # Set telemetry bus for each component
        self._set_telemetry_bus()

        # Initialize numpy array associated to each component
        self._vectorize_components()

        # Partition vector_regions into approximately equally sized batches
        if self.enable_parallel:
            self._parallelize_vector_regions()
        else:
            self.num_jobs = 1

        # Assign prngs to each vectorized region and set random seeds
        self._seed_regions()

        # Initialise components, such as disease model, movement model, interventions etc
        self._initial_conditions(offset)

    def run(self):
        """Run the simulation"""

        ticks_in_week = self.clock.ticks_in_week
        ticks_in_day = self.clock.ticks_in_day
        offset = self.clock.ticks_through_week()

        # Notify telemetry bus of simulation data
        self._output_data_intial()

        log.info("Simulating outbreak on " + str(self.num_jobs) + " CPU(s)...")

        # Start the main loop
        for day in self.clock:

            self.telemetry_bus.publish("sim.time", self.clock)

            self.travel_model.dynamics(self, day, ticks_in_day,
                                                self.health_model.facemask_transmission_multiplier,
                                                self.health_model.mutation_matrix,
                                                self.enable_parallel,
                                                self.num_jobs,
                                                self.vector_region_batches)

            if self.enable_parallel:
                Parallel(n_jobs=self.num_jobs, backend="threading",
                         verbose=0)(delayed(self.simulate_day)(vector_region_batch, day, offset,
                                                               ticks_in_day, ticks_in_week)
                                    for vector_region_batch in self.vector_region_batches)
            else:
                self.simulate_day(self.vector_regions, day, offset, ticks_in_day, ticks_in_week)

            self._output_data_update(day)

        # Notify the telemetry bus that the simulation has ended
        self.telemetry_bus.publish("simulation.end")

    def _update(self, vector_region, t):
        """Update the agents."""

        # Update health
        self.health_model.update(vector_region, t)

        # Update current location and facemasks
        self.movement_model.update(vector_region, t)

    def _output_data_intial(self):
        """Published intial strain count data to the telemetry bus"""

        self.day_to_date = {}

        self.cumulative_deaths_old = {vr.name: 0 for vr in self.vector_regions}
        self.total_deaths = 0
        self.daily_deaths = defaultdict(dict)
        self.cumulative_deaths = defaultdict(dict)
        self.reporter_deaths = defaultdict(dict)

        if self.config['reporters'] is not None:

            population_sizes =\
                [vector_region.number_of_agents for vector_region in self.vector_regions]
            region_names =\
                [vector_region.name for vector_region in self.vector_regions]
            region_names_to_super_regions =\
                {vr.name: vr.super_region for vr in self.vector_regions}

            self.telemetry_bus.publish("strain_counts.initialize", self.number_of_regions,
                                                                   self.number_of_strains,
                                                                   region_names,
                                                                   population_sizes)

            regions_omitted = self.config['regions_omitted_from_death_counts']

            self.telemetry_bus.publish("deaths.initialize", region_names_to_super_regions,
                                                            regions_omitted)

            self.telemetry_bus.publish("vector_world.data", self.vector_world,
                                                            self.health_model.number_of_strains)

            region_to_omit = self.config['regions_omitted_from_death_counts']

            self.telemetry_bus.publish("death_counts.initialize", self.number_of_regions,
                                                                  self.number_of_strains,
                                                                  region_names,
                                                                  population_sizes,
                                                                  region_to_omit)

    def _output_data_update(self, day):
        """Published strain count updates to the telemetry bus"""

        date = self.clock.iso8601()
        self.day_to_date[day] = date

        total_deaths = 0
        for vector_region in self.vector_regions:
            if vector_region.name not in self.config['regions_omitted_from_death_counts']:
                cumulative_deaths_new = 0
                cumulative_deaths_new +=\
                    self.count_dead(
                        c_int64(vector_region.number_of_agents),
                        c_void_p(vector_region.current_disease.ctypes.data),
                    )
                total_deaths += cumulative_deaths_new
                self.daily_deaths[date][vector_region.name] =\
                    int((1 / self.vector_world.scale_factor) *\
                    (cumulative_deaths_new - self.cumulative_deaths_old[vector_region.name]))
                self.cumulative_deaths[date][vector_region.name] =\
                    int((1 / self.vector_world.scale_factor) * cumulative_deaths_new)
                self.cumulative_deaths_old[vector_region.name] = cumulative_deaths_new

        self.total_deaths = int((1 / self.vector_world.scale_factor) * total_deaths)

        if ("pygame_coords.PygameCoords" in self.config['reporters']):

            current_location =\
                {vr.id: vr.current_location for vr in self.vector_regions}
            current_infectiousness =\
                {vr.id: vr.current_infectiousness for vr in self.vector_regions}

            self.telemetry_bus.publish("data.update", self.clock, current_location,
                                       current_infectiousness)

        if ("plot.PlotInfected" in self.config['reporters']) or\
            ("csv.StrainCounts" in self.config['reporters']) or\
            ("pygame_shapes.PygameShapes" in self.config['reporters']):

            infections =\
                np.zeros((self.number_of_regions * self.number_of_strains), dtype=np.int64)

            for vector_region in self.vector_regions:
                self.collect_telemetry_data(
                    c_int64(vector_region.number_of_agents),
                    c_int64(self.number_of_strains),
                    c_int64(vector_region.id),
                    c_void_p(vector_region.current_strain.ctypes.data),
                    c_void_p(infections.ctypes.data)
                )
            infections = ((1 / self.vector_world.scale_factor) * infections).astype(np.int64)
            self.telemetry_bus.publish("strain_counts.update", self.clock, infections)

        if ("csv.DeathCounts" in self.config['reporters']):

            self.telemetry_bus.publish("death_counts.update", self.clock, self.daily_deaths[date])

        if ("plot.PlotDeaths" in self.config['reporters']):

            self.telemetry_bus.publish("deaths.update", self.clock, self.daily_deaths[date],
                                        self.cumulative_deaths[date])

    def calculate_cost(self, policy):
        """Calculates the final cost of the pandemic, to be used for policy optimization"""

        log.info("Total deaths: %d", self.total_deaths)

        return self.total_deaths

    def calculate_error(self):
        """Calculates rolling 21 day average of daily deaths and historical deaths and calculates
        difference"""

        omitted_region_names = []
        for vector_region in self.vector_regions:
            if (vector_region.super_region in\
                self.config['super_regions_omitted_from_deaths_counts']) or\
                (vector_region.name in\
                self.config['regions_omitted_from_death_counts']):
                omitted_region_names.append(vector_region.name)
        
        average_deaths_dict_simulation = {}

        # Calculate 21-day rolling average of daily deaths
        total_num_days = self.clock.simulation_length_days
        for day in range(total_num_days):
            min_day = max(0, day - 10)
            max_day = min(total_num_days, day + 10)
            num_days = max_day - min_day
            average_deaths = 0
            for day_sample in range(min_day, max_day):
                for vector_region in self.vector_regions:
                    if vector_region.name not in omitted_region_names:
                        average_deaths +=\
                            self.daily_deaths[self.day_to_date[day_sample]][vector_region.name]
            average_deaths /= num_days
            average_deaths_dict_simulation[self.day_to_date[day]] = average_deaths

        included_region_names = []
        for vector_region in self.vector_regions:
            if vector_region.name not in omitted_region_names:
                included_region_names.append(vector_region.name)

        # Load historical data
        if self.average_deaths_dict_historical is None:
            historical_daily_deaths =\
                {date: average_deaths_dict_simulation[date] for\
                date in average_deaths_dict_simulation}
            if self.historical_data_filepath is not None:
                with open(self.historical_data_filepath, newline='') as csvfile:
                    next(csvfile)
                    deaths_data = csv.reader(csvfile, delimiter=',')
                    for row in deaths_data:
                        iso2 = str(row[1])
                        if iso2 in included_region_names:
                            date = datetime.strptime(str(row[0]), '%Y-%m-%d').strftime('%d/%m/%Y%Z')
                            if date in historical_daily_deaths:
                                historical_daily_deaths[date] += int(row[6])
                            else:
                                historical_daily_deaths[date] = int(row[6])

            self.average_deaths_dict_historical = {}

            # Calculate 21-day rolling average
            total_num_days = self.clock.simulation_length_days
            for day in range(total_num_days):
                min_day = max(0, day - 10)
                max_day = min(total_num_days, day + 10)
                num_days = max_day - min_day
                average_deaths = 0
                for day_sample in range(min_day, max_day):
                    average_deaths += historical_daily_deaths[self.day_to_date[day_sample]]
                average_deaths /= num_days
                self.average_deaths_dict_historical[self.day_to_date[day]] = average_deaths

        # Calculate L^2 difference
        for date in average_deaths_dict_simulation:
            self.error += (average_deaths_dict_simulation[date] -\
                          self.average_deaths_dict_historical[date])**2
        self.error = np.sqrt(self.error)

        log.info("Simulation Error: %f", self.error)
