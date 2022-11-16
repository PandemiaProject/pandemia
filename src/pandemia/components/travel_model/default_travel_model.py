"""Default regional mixing model"""

import logging
import numpy as np
import copy
from joblib import Parallel, delayed

from ctypes import c_void_p, c_double, c_int, cdll

from ..travel_model import TravelModel

log = logging.getLogger("default_travel_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultTravelModel(TravelModel):
    """Default model of mixing between regions. Each day, a number of agents are selected from each
    region to travel to each other region. The numbers travelling between regions are given by the
    matrix agents_travelling_matrix. The destination of the travellers is recorded using the array
    current_region. Agents outside their home region are typically ignored by other components
    while they are travelling. Only uninfected agents can travel, but such agents can become
    infected whlie travelling. The probability that they become infected is given in terms of the
    average infectiousness of their destination region. The border closure intervention is also
    implemented here."""

    def __init__(self, config, scale_factor, number_of_strains,
                 vector_world):
        """Initialize component"""
        super().__init__(config, scale_factor)

        lib = cdll.LoadLibrary("./src/pandemia/components/travel_model/"
                                "default_travel_model_functions.dll")

        self.transmission_out = lib.transmission_out
        self.transmission_in = lib.transmission_in
        self.close_borders = lib.close_borders
        self.determine_travellers = lib.determine_travellers

        self.transmission_out.restype = None
        self.transmission_in.restype = None
        self.close_borders.restype = None
        self.determine_travellers.restype = None

        self.scale_factor      = scale_factor
        self.number_of_regions = vector_world.number_of_regions
        self.number_of_strains = number_of_strains

        self.travel_multiplier = config['travel_transmission_multiplier']
        self.interpolation     = config['interpolation']

        self.baseline_agents_travelling_matrix = vector_world.travel_matrix
        for r in range(self.number_of_regions):
            assert self.baseline_agents_travelling_matrix[r][r] == 0

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        vector_region.current_border_closure_multiplier = 1.0
        vector_region.current_region = np.full((vector_region.number_of_agents),
                                                vector_region.id, dtype=int)

    def initial_conditions(self, sim):
        """Initial regional mixing"""

        self.beta = sim.health_model.beta

        self.baseline_agents_travelling_matrix =\
            self.interpolate_matrix(self.baseline_agents_travelling_matrix, sim.vector_world,
                                    self.interpolation)

        # Flatten arrays
        self.baseline_agents_travelling_matrix =\
            self.baseline_agents_travelling_matrix.flatten()

    def _out(self, vector_region_batch, agents_travelling_matrix, sum_f_by_strain,
             transmission_force, facemask_transmission_multiplier):
        """Update record of who is travelling abroad for a batch of regions"""

        for vector_region in vector_region_batch:
            self.determine_travellers(
                c_int(vector_region.number_of_agents),
                c_int(vector_region.id),
                c_int(self.number_of_regions),
                c_void_p(vector_region.current_disease.ctypes.data),
                c_void_p(vector_region.current_strain.ctypes.data),
                c_void_p(vector_region.current_region.ctypes.data),
                c_void_p(agents_travelling_matrix.ctypes.data),
                c_void_p(vector_region.random_state.ctypes.data)
            )
            self.transmission_out(
                c_int(vector_region.number_of_agents),
                c_int(self.number_of_strains),
                c_int(vector_region.id),
                c_void_p(self.beta.ctypes.data),
                c_double(facemask_transmission_multiplier),
                c_double(self.travel_multiplier),
                c_double(vector_region.current_region_transmission_multiplier),
                c_void_p(vector_region.current_region.ctypes.data),
                c_void_p(vector_region.current_infectiousness.ctypes.data),
                c_void_p(vector_region.current_strain.ctypes.data),
                c_void_p(vector_region.current_facemask.ctypes.data),
                c_void_p(sum_f_by_strain.ctypes.data),
                c_void_p(transmission_force.ctypes.data),
                c_void_p(vector_region.random_state.ctypes.data)
            )

    def _in(self, sim, vector_region_batch, sum_f_by_strain, transmission_force, mutation_matrix,
            facemask_transmission_multiplier, day, ticks_in_day):
        """Calculate who gets infected for a batch of regions"""

        for vector_region in vector_region_batch:
            self.transmission_in(
                c_int(self.number_of_regions),
                c_int(self.number_of_strains),
                c_int(vector_region.number_of_agents),
                c_int(vector_region.id),
                c_void_p(vector_region.current_facemask.ctypes.data),
                c_void_p(vector_region.current_region.ctypes.data),
                c_double(facemask_transmission_multiplier),
                c_void_p(sum_f_by_strain.ctypes.data),
                c_void_p(vector_region.current_sigma_immunity_failure.ctypes.data),
                c_void_p(vector_region.infection_event.ctypes.data),
                c_void_p(transmission_force.ctypes.data),
                c_void_p(mutation_matrix.ctypes.data),
                c_void_p(vector_region.random_state.ctypes.data)
            )
            time = (day * ticks_in_day) +\
                    vector_region.prng.random_randrange(ticks_in_day)
            sim.health_model.infect_wrapper(vector_region, time)
            sim.health_model.update(vector_region, day * ticks_in_day)

    def dynamics(self, sim, day, ticks_in_day,
                 facemask_transmission_multiplier,
                 mutation_matrix, enable_parallel,
                 num_jobs, vector_region_batches):
        """Changes to regional mixing model"""

        agents_travelling_matrix =\
            np.zeros((self.number_of_regions, self.number_of_regions), dtype=int)

        # Reset record of who is travelling abroad and apply border closure if necessary
        for vector_region in sim.vector_regions:
            vector_region.current_region = np.full((vector_region.number_of_agents),
                                                   vector_region.id, dtype=int)
            self.close_borders(
                c_int(self.number_of_regions),
                c_int(vector_region.id),
                c_double(self.scale_factor),
                c_double(vector_region.current_border_closure_multiplier),
                c_void_p(agents_travelling_matrix.ctypes.data),
                c_void_p(self.baseline_agents_travelling_matrix.ctypes.data)
            )

        # Update record of who is travelling abroad and calculate transmission force for each region
        sum_f_by_strain = np.zeros((self.number_of_regions * self.number_of_strains), dtype=float)
        transmission_force = np.ones((self.number_of_regions), dtype=float)
        if enable_parallel:
            Parallel(n_jobs=num_jobs, backend="threading",
                     verbose=0)(delayed(self._out)(vector_region_batch, agents_travelling_matrix,
                     sum_f_by_strain, transmission_force, facemask_transmission_multiplier)
                     for vector_region_batch in vector_region_batches)
        else:
            self._out(sim.vector_regions, agents_travelling_matrix,
                      sum_f_by_strain, transmission_force, facemask_transmission_multiplier)

        # Calculate who gets infected and infect them
        if enable_parallel:
            Parallel(n_jobs=num_jobs, backend="threading",
                     verbose=0)(delayed(self._in)(sim, vector_region_batch,
                     sum_f_by_strain, transmission_force, mutation_matrix,
                     facemask_transmission_multiplier, day, ticks_in_day)
                     for vector_region_batch in vector_region_batches)
        else:
            self._in(sim, sim.vector_regions, sum_f_by_strain, transmission_force, mutation_matrix,
                     facemask_transmission_multiplier, day, ticks_in_day)

    def interpolate_matrix(self, baseline_agents_travelling_matrix, vector_world, interpolation):
        """Artificially inflate number of travellers for testing purposes or otherwise"""

        ids_to_population_sizes =\
            {vr.id: int(vr.number_of_agents * (1 / self.scale_factor))\
                        for vr in vector_world.vector_regions}

        # Adjust matrix using interpolation parameter
        row_sums = baseline_agents_travelling_matrix.sum(axis=1)
        for i in range(self.number_of_regions):
            for j in range(self.number_of_regions):
                if baseline_agents_travelling_matrix[i][j] > 0:
                    max_proportion_travelling =\
                        baseline_agents_travelling_matrix[i][j] / row_sums[i]
                    proportion_travelling =\
                        baseline_agents_travelling_matrix[i][j] / ids_to_population_sizes[i]
                    interpolated_proportion =\
                        (proportion_travelling * (1 - interpolation)) +\
                            (max_proportion_travelling * interpolation)
                    baseline_agents_travelling_matrix[i][j] =\
                        (ids_to_population_sizes[i] * interpolated_proportion)

        return baseline_agents_travelling_matrix
