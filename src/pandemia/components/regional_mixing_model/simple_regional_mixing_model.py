"""Simple regional mixing model"""

import logging
import numpy as np
import copy
from joblib import Parallel, delayed

from ctypes import c_void_p, pointer, c_double, c_int, cdll

from pandemia.components.regional_mixing_model import RegionalMixingModel

log = logging.getLogger("simple_regional_mixing_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class SimpleRegionalMixingModel(RegionalMixingModel):
    """Simple model of regional mixing"""

    def __init__(self, config, scale_factor, number_of_strains,
                 vector_world, enable_ctypes):
        """Initialize component"""
        super().__init__(config, scale_factor)

        self.enable_ctypes = enable_ctypes

        assert self.enable_ctypes, "Must use ctypes, Python functions out-of-date"

        if self.enable_ctypes:

            lib = cdll.LoadLibrary("./src/pandemia/components/regional_mixing_model/"
                                   "simple_regional_mixing_model_functions.dll")

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

        # If using ctypes, flatten arrays accordingly
        if self.enable_ctypes:
            self.baseline_agents_travelling_matrix =\
                self.baseline_agents_travelling_matrix.flatten()

    def dynamics(self, sim, day, ticks_in_day,
                 facemask_transmission_multiplier,
                 mutation_matrix, enable_parallel,
                 num_jobs, vector_region_batches):
        """Changes agent locations"""

        if self.enable_ctypes:
            self.dynamics_c(sim, day, ticks_in_day,
                            facemask_transmission_multiplier,
                            mutation_matrix, enable_parallel,
                            num_jobs, vector_region_batches)
        else:
            self.dynamics_python(sim, day, ticks_in_day,
                                 facemask_transmission_multiplier,
                                 mutation_matrix)

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
                c_void_p(vector_region.random_state.ctypes.data),
                pointer(vector_region.random_p)
            )
            self.transmission_out(
                c_int(vector_region.number_of_agents),
                c_int(self.number_of_strains),
                c_int(vector_region.id),
                c_double(self.beta),
                c_double(facemask_transmission_multiplier),
                c_double(self.travel_multiplier),
                c_double(vector_region.current_region_transmission_multiplier),
                c_void_p(vector_region.current_region.ctypes.data),
                c_void_p(vector_region.current_infectiousness.ctypes.data),
                c_void_p(vector_region.current_strain.ctypes.data),
                c_void_p(vector_region.current_facemask.ctypes.data),
                c_void_p(sum_f_by_strain.ctypes.data),
                c_void_p(transmission_force.ctypes.data),
                c_void_p(vector_region.random_state.ctypes.data),
                pointer(vector_region.random_p)
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
                c_void_p(vector_region.random_state.ctypes.data),
                pointer(vector_region.random_p)
            )
            time = (day * ticks_in_day) +\
                    vector_region.prng.random_randrange(ticks_in_day)
            sim.health_model.infect_c(vector_region, time)
            sim.health_model.update_c(vector_region, day * ticks_in_day)

    def dynamics_c(self, sim, day, ticks_in_day,
                   facemask_transmission_multiplier,
                   mutation_matrix, enable_parallel,
                   num_jobs, vector_region_batches):
        """Changes to regional mixing"""

        agents_travelling_matrix = copy.deepcopy(self.baseline_agents_travelling_matrix)

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

    def dynamics_python(self, sim, day, ticks_in_day,
                        facemask_transmission_multiplier,
                        mutation_matrix):
        """Changes to regional mixing"""

        # Reset record of who is travelling abroad
        for vector_region in sim.vector_regions:
            vector_region.current_region = np.full((vector_region.number_of_agents),
                                                   vector_region.id, dtype=int)

        # Apply border closure intervention if necessary
        agents_travelling_matrix = copy.deepcopy(self.baseline_agents_travelling_matrix)
        route_suppressed = np.zeros((self.number_of_regions, self.number_of_regions), dtype=int)
        for vector_region in sim.vector_regions:
            if vector_region.border_closure_intervention == 1:
                for other_vector_region in sim.vector_regions:
                    r1 = vector_region.id
                    r2 = other_vector_region.id
                    if r1 != r2:
                        route_suppressed[r1][r2] = 1
        for r1 in range(self.number_of_regions):
            for r2 in range(self.number_of_regions):
                if route_suppressed[r1][r2] == 1:
                    agents_travelling_matrix[r1][r2] =\
                        int(self.baseline_agents_travelling_matrix[r1][r2] *\
                            self.border_closure_factor)

        # Update record of who is travelling abroad
        for vector_region in sim.vector_regions:
            r1 = vector_region.id
            total_num_to_travel = sum(agents_travelling_matrix[r1])
            if total_num_to_travel > 0:
                eligible_to_travel = [n for n in range(vector_region.number_of_agents)
                                      if (vector_region.current_strain[n] == -1 and
                                      vector_region.current_disease[n] < 1.0)]
                for r2 in range(self.number_of_regions):
                    num_to_travel = agents_travelling_matrix[r1][r2]
                    if num_to_travel > 0:
                        num = min(num_to_travel, len(eligible_to_travel))
                        travellers = vector_region.prng.random_sample(eligible_to_travel, num)
                        for n in travellers:
                            eligible_to_travel.remove(n)
                            vector_region.current_region[n] = r2

        # Calculate average transmission force for each region
        transmission_force = np.zeros((self.number_of_regions, self.number_of_regions), dtype=float)
        sum_f_by_strain = np.zeros((self.number_of_regions, self.number_of_strains), dtype=float)
        average_f = np.zeros((self.number_of_regions), dtype=float)
        for vector_region in sim.vector_regions:
            for n in range(vector_region.number_of_agents):
                if vector_region.current_region[n] == vector_region.id:
                    if vector_region.current_strain[n] != -1:
                        facemask_multiplier = 1 + vector_region.current_facemask[n] *\
                                              (facemask_transmission_multiplier - 1)
                        region_multiplier = vector_region.current_region_transmission_multiplier
                        infectiousness = vector_region.current_infectiousness[n]
                        f = region_multiplier * facemask_multiplier * infectiousness
                        sum_f_by_strain[vector_region.id][vector_region.current_strain[n]] += f
            average_f[vector_region.id] = sum(sum_f_by_strain[vector_region.id])/\
                                              vector_region.number_of_agents
        for r1 in range(self.number_of_regions):
            for r2 in range(self.number_of_regions):
                num_contacted = self.contacts_matrix[r1][r2]
                travel_multiplier = self.travel_transmission_multiplier[r1][r2]
                duration = self.contact_ticks_matrix[r1][r2]
                transmission_force[r1][r2] = 1 - ((1 - travel_multiplier * average_f[r2]) **\
                                             (duration * num_contacted))

        # Calculate who gets infected
        for vector_region in sim.vector_regions:
            r1 = vector_region.id
            for n in range(vector_region.number_of_agents):
                if vector_region.current_region != r1:
                    r2 = vector_region.current_region.id
                    facemask_multiplier = 1 + vector_region.current_facemask[n] *\
                                              (facemask_transmission_multiplier - 1)
                    if vector_region.prng.random_float(1.0) <\
                        facemask_multiplier * transmission_force[r1][r2]:
                        strains = [s for s in range(self.number_of_strains)]
                        weights = sum_f_by_strain[r2]
                        s1 = vector_region.prng.random_choices(strains, weights, 1)[0]
                        if vector_region.prng.random_float(1.0) <\
                            vector_region.current_sigma_immunity_failure[n][s1]:
                            vector_region.infection_event[n] =\
                                self._mutate(s1, mutation_matrix, vector_region.prng)

        # Infect
        for vector_region in sim.vector_regions:
            time = (day * ticks_in_day) +\
                    vector_region.prng.random_randrange(ticks_in_day)
            sim.health_model.infect_python(vector_region, time)
            sim.health_model.update_python(vector_region, day * ticks_in_day)

    def _mutate(self, s1, mutation_matrix, prng):
        """Mutates strain"""

        weights = mutation_matrix[s1]
        strains = [s2 for s2 in range(len(weights))]
        return prng.random_choices(strains, weights, 1)[0]
