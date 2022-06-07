"""Simulates an epidemic"""

import logging
import numpy as np
import os
from datetime import datetime
from ctypes import c_int, c_void_p, cdll
from joblib import Parallel, delayed

from pandemia.random_tools import Random

import uuid

from pandemia.version import VERSION

log = logging.getLogger('sim')

# TODO
# - define environment and dependencies etc
# - generate documentation
# - compile for different operating systems and test

# TESTING
# - tests and assertions, e.g. assertions for all presets and vaccines, e.g. for
#   number_of_rho_immunity_outcomes (assert that number_of_rho_immunity_outcomes is the same for
#   both health model and vaccination model etc...)

# SCENARIOS
# - configure other scenarios (i.e. world objects) using data, e.g. the Luxembourg data...

# PERFORMANCE
# - need to check don't overflow stack when allocating memory in the C functions, if array length
#   (e.g. population size) too big etc...? (move those arrays outside of stack?)
# - optimize runtime and memory usage (padded arrays can be compressed with offsets, for example)

# FEATURES
# - variable resolution, meaning a different scale factor for each region(?)
# - rigourous statistical validation methods?
# - include a regional/location transmission factor intervention respresenting social distancing?
# - impact of hospitalization and capacity?
# - testing and contact tracing delays?
# - improved immunity subsampling approach?
# - fast forward remaining sim if total infected is zero?

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
                 regional_mixing_model,
                 input_model,
                 telemetry_bus):

        # Static info
        self.pandemia_version = VERSION
        self.created_at       = datetime.now()
        self.run_id           = uuid.uuid4().hex

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
        self.regional_mixing_model = regional_mixing_model
        self.input_model = input_model

        self.number_of_strains = self.health_model.number_of_strains
        self.number_of_vaccines = self.vaccination_model.number_of_vaccines
        self.number_of_regions = vector_world.number_of_regions

        # Enable ctypes
        self.enable_ctypes = self.config['enable_ctypes']

        if self.enable_ctypes:

            lib = cdll.LoadLibrary("./src/pandemia/simulator_functions.dll")

            self.collect_telemetry_data = lib.collect_telemetry_data

            self.collect_telemetry_data.restype = None

        # Parallel processing
        self.enable_parallel       = self.config['enable_parallel']
        self.num_jobs              = self.config['num_jobs']
        self.vector_region_batches = None

    def _set_telemetry_bus(self):
        """Assigns telemetry bus to each component"""

        self.seasonal_effects_model.set_telemetry_bus(self.telemetry_bus)
        self.health_model.set_telemetry_bus(self.telemetry_bus)
        self.movement_model.set_telemetry_bus(self.telemetry_bus)
        self.hospitalization_and_death_model.set_telemetry_bus(self.telemetry_bus)
        self.testing_and_contact_tracing_model.set_telemetry_bus(self.telemetry_bus)
        self.vaccination_model.set_telemetry_bus(self.telemetry_bus)
        self.regional_mixing_model.set_telemetry_bus(self.telemetry_bus)
        self.input_model.set_telemetry_bus(self.telemetry_bus)

    def _vectorize_components(self):
        """Initialize the numpy arrays associated to each component"""

        for vector_region in self.vector_regions:

            self.regional_mixing_model.vectorize_component(vector_region)
            self.input_model.vectorize_component(vector_region)
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
            self.num_jobs = min(max(number_of_cpus + 1 + self.num_jobs, 1), number_of_cpus)
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
            # Each vector region also gets a state and counter for the prng used inside C
            np.random.seed(seed)
            vector_region.random_state =\
                np.random.randint(0, ULONG_MAX + 1, size=16, dtype=np.uint64)
            vector_region.random_p = c_int(0)

    def _initial_conditions(self, offset):
        """Initialize the various submodels"""

        self.regional_mixing_model.initial_conditions(self)

        for vector_region in self.vector_regions:

            self.input_model.initial_conditions(vector_region)
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

            self.input_model.dynamics(vector_region, day)
            self.seasonal_effects_model.dynamics(vector_region, day)

            for tick in range(ticks_in_day):

                t = (ticks_in_day * day) + tick

                self.health_model.dynamics(vector_region, t)
                self.movement_model.dynamics(vector_region, t, offset, ticks_in_week)
                self.hospitalization_and_death_model.dynamics(vector_region)

                self._update(vector_region, t)

            self.testing_and_contact_tracing_model.dynamics(vector_region, day)
            self.vaccination_model.dynamics(vector_region, day, ticks_in_day)

    def run(self):
        """Run the simulation"""

        # Set the correct time
        self.clock.reset()
        ticks_in_week = self.clock.ticks_in_week
        ticks_in_day = self.clock.ticks_in_day
        offset = self.clock.ticks_through_week()

        # Set telemetry bus for each component
        self._set_telemetry_bus()

        # Initialize numpy array associated to each component
        self._vectorize_components()

        # Partition vector_regions into approximately equally sized batches
        if self.enable_parallel and self.enable_ctypes:
            self._parallelize_vector_regions()

        # Assign prngs to each vectorized region and set random seeds
        self._seed_regions()

        # Initialise components, such as disease model, movement model, interventions etc
        self._initial_conditions(offset)

        # Notify telemetry bus of simulation data
        self._publish_telemetry_data_intial()

        log.info("Simulating outbreak on " + str(self.num_jobs) + " CPUs...")

        # Start the main loop
        for day in self.clock:

            self.telemetry_bus.publish("sim.time", self.clock)

            self.regional_mixing_model.dynamics(self, day, ticks_in_day,
                                                self.health_model.facemask_transmission_multiplier,
                                                self.health_model.mutation_matrix,
                                                self.enable_parallel,
                                                self.num_jobs,
                                                self.vector_region_batches)

            if self.enable_parallel and self.enable_ctypes:
                Parallel(n_jobs=self.num_jobs, backend="threading",
                         verbose=0)(delayed(self.simulate_day)(vector_region_batch, day, offset,
                                                               ticks_in_day, ticks_in_week)
                                    for vector_region_batch in self.vector_region_batches)
            else:
                self.simulate_day(self.vector_regions, day, offset, ticks_in_day, ticks_in_week)

            self._publish_telemetry_data_update()

        # Notify the telemetry bus that the simulation has ended
        self.telemetry_bus.publish("simulation.end")

    def _update(self, vector_region, t):
        """Update the agents."""

        # Update health
        self.health_model.update(vector_region, t)

        # Update current location and facemasks
        self.movement_model.update(vector_region, t)

    def _publish_telemetry_data_intial(self):
        """Published intial strain count data to the telemetry bus"""

        population_sizes = [vector_region.number_of_agents for vector_region in self.vector_regions]
        region_names = [vector_region.name for vector_region in self.vector_regions]

        self.telemetry_bus.publish("strain_counts.initialize", self.number_of_regions,
                                                               self.number_of_strains,
                                                               region_names,
                                                               population_sizes)

        self.telemetry_bus.publish("vector_world.data", self.vector_world,
                                                        self.health_model.number_of_strains)

    def _publish_telemetry_data_update(self):
        """Published strain count updates to the telemetry bus"""

        strain_counts = np.zeros((self.number_of_regions * self.number_of_strains), dtype=int)
        for vector_region in self.vector_regions:
            if self.enable_ctypes:
                self.collect_telemetry_data(
                    c_int(vector_region.number_of_agents),
                    c_int(self.number_of_strains),
                    c_int(vector_region.id),
                    c_void_p(vector_region.current_strain.ctypes.data),
                    c_void_p(strain_counts.ctypes.data)
                )
            else:
                for n in range(vector_region.number_of_agents):
                    if vector_region.current_strain[n] != -1:
                        strain_counts[(vector_region.id * self.number_of_strains)
                                      + vector_region.current_strain[n]] += 1
        strain_counts = int(1 / self.vector_world.scale_factor) * strain_counts
        self.telemetry_bus.publish("strain_counts.update", self.clock, strain_counts)

    def total_deaths(self):
        """Calculates final death count"""

        total_deaths = 0
        for vector_region in self.vector_regions:
            for n in range(vector_region.number_of_agents):
                if vector_region.current_disease[n] == 1.0:
                    total_deaths += 1
        scale_factor = self.vector_world.scale_factor
        total_deaths = int((1 / scale_factor) * total_deaths)

        return total_deaths
