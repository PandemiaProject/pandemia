"""Safir regional mixing model"""

import logging
import numpy as np
import copy
from joblib import Parallel, delayed

from ctypes import c_void_p, c_double, c_int, cdll

from ..travel_model import TravelModel
log = logging.getLogger("safir_travel_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class SafirTravelModel(TravelModel):
    """Version of the default model of mixing between regions adapted to the safir health model.

    Each day, a number of agents are selected from each region to travel to each other region. The
    numbers travelling between regions are given by the matrix agents_travelling_matrix. The
    destination of the travellers is recorded using the array current_region. Agents outside their
    home region are typically ignored by other components while they are travelling. Only uninfected
    agents can travel, but such agents can become infected whlie travelling. The probability that
    they become infected is given in terms of the average infectiousness of their destination
    region. The border closure intervention is also implemented here.

    Parameters:
    -----------
    config : Config
        A Pandemia Config object. A sub-config of the full config, containing the data used to
        configure this component.
    scale_factor : float
        The scale factor, coming from the full config.
    number_of_strains : int
        The number of strains appearing the model.
    vector_world : VectorWorld
        A Pandemia VectorWorld object, containing the relevant regions.
    """

    def __init__(self, config, scale_factor, number_of_strains,
                 vector_world):
        """Initialize component"""
        super().__init__(config, scale_factor)

        self.transmission_out = self.lib.safir_3_transmission_out
        self.transmission_in = self.lib.safir3_transmission_in
        self.close_borders = self.lib.close_borders
        self.determine_travellers = self.lib.determine_travellers

        self.transmission_out.restype = None
        self.transmission_in.restype = None
        self.close_borders.restype = None
        self.determine_travellers.restype = None

        self.scale_factor      = scale_factor
        self.number_of_regions = vector_world.number_of_regions
        self.number_of_strains = number_of_strains

        self.travel_multiplier = config['travel_transmission_multiplier']
        self.interpolation     = config['interpolation']
        self.beta              = np.array(config['beta'], dtype=np.float64)

        self.baseline_agents_travelling_matrix = vector_world.travel_matrix
        for r in range(self.number_of_regions):
            assert self.baseline_agents_travelling_matrix[r][r] == 0

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        vector_region.current_border_closure_multiplier = 1.0
        vector_region.current_region = np.full((vector_region.number_of_agents),
                                                vector_region.id, dtype=np.int32)

    def initial_conditions(self, sim):
        """Initial regional mixing"""

        self.baseline_agents_travelling_matrix =\
            self._interpolate_matrix(self.baseline_agents_travelling_matrix, sim.vector_world,
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
                c_int(vector_region.id),
                c_void_p(self.beta.ctypes.data),
                c_double(facemask_transmission_multiplier),
                c_double(self.travel_multiplier),
                c_double(vector_region.current_region_transmission_multiplier),
                c_void_p(vector_region.current_region.ctypes.data),
                c_void_p(vector_region.current_infectiousness.ctypes.data),
                c_void_p(vector_region.ef_transmission.ctypes.data),
                c_void_p(vector_region.current_strain.ctypes.data),
                c_void_p(vector_region.current_facemask.ctypes.data),
                c_void_p(transmission_force.ctypes.data)
            )

    def _in(self, sim, vector_region_batch, sum_f_by_strain, transmission_force, mutation_matrix,
            facemask_transmission_multiplier, day, ticks_in_day):
        """Calculate who gets infected for a batch of regions"""

        for vector_region in vector_region_batch:
            self.transmission_in(
                c_int(vector_region.number_of_agents),
                c_int(vector_region.id),
                c_void_p(vector_region.current_facemask.ctypes.data),
                c_void_p(vector_region.current_region.ctypes.data),
                c_double(facemask_transmission_multiplier),
                c_void_p(vector_region.ef_infection.ctypes.data),
                c_void_p(vector_region.infection_event.ctypes.data),
                c_void_p(transmission_force.ctypes.data),
                c_void_p(vector_region.random_state.ctypes.data)
            )

            time = (day * ticks_in_day) +\
                    vector_region.prng.random_randrange(ticks_in_day)

            indexes = np.argwhere(vector_region.infection_event == 0).flatten()
            vector_region.t_for_recovery[indexes] = time + sim.health_model.infectious_period_ticks
            vector_region.start_t[indexes] = time
            vector_region.current_strain[indexes] = 0
            vector_region.infection_event =\
                np.full(vector_region.number_of_agents, -1, dtype=np.int32)

    def dynamics(self, sim, day, ticks_in_day,
                 facemask_transmission_multiplier,
                 mutation_matrix, enable_parallel,
                 num_jobs, vector_region_batches):
        """Changes to regional mixing model"""

        agents_travelling_matrix =\
            np.zeros((self.number_of_regions, self.number_of_regions), dtype=np.int32)

        # Reset record of who is travelling abroad and apply border closure if necessary
        for vector_region in sim.vector_regions:
            vector_region.current_region = np.full((vector_region.number_of_agents),
                                                   vector_region.id, dtype=np.int32)
            self.close_borders(
                c_int(self.number_of_regions),
                c_int(vector_region.id),
                c_double(self.scale_factor),
                c_double(vector_region.current_border_closure_multiplier),
                c_void_p(agents_travelling_matrix.ctypes.data),
                c_void_p(self.baseline_agents_travelling_matrix.ctypes.data)
            )

        # Update record of who is travelling abroad and calculate transmission force for each region
        sum_f_by_strain =\
            np.zeros((self.number_of_regions * self.number_of_strains), dtype=np.float64)
        transmission_force =\
            np.ones((self.number_of_regions), dtype=np.float64)
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

    def _interpolate_matrix(self, baseline_agents_travelling_matrix, vector_world, interpolation):
        """Artificially inflate number of travellers for testing purposes or otherwise"""

        ids_to_population_sizes =\
            {vr.id: np.int32(vr.number_of_agents * (1 / self.scale_factor))\
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
