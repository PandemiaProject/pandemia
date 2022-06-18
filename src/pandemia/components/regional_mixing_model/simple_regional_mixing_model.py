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

        if self.enable_ctypes:

            lib = cdll.LoadLibrary("./src/pandemia/components/regional_mixing_model/"
                                   "simple_regional_mixing_model_functions.dll")

            self.calculate_average_f = lib.calculate_average_f
            self.average_transmission_force = lib.average_transmission_force
            self.travel = lib.travel
            self.suppress_routes = lib.suppress_routes
            self.close_borders = lib.close_borders
            self.determine_travellers = lib.determine_travellers

            self.calculate_average_f.restype = None
            self.average_transmission_force.restype = None
            self.travel.restype = None
            self.suppress_routes.restype = None
            self.close_borders.restype = None
            self.determine_travellers.restype = None

        self.scale_factor      = scale_factor
        self.number_of_regions = vector_world.number_of_regions
        self.number_of_strains = number_of_strains

        self.border_closure_factor          = config['border_closure_factor']
        self.travel_transmission_multiplier = self.config['travel_transmission_multiplier']

        self.travel_transmission_multiplier =\
            np.full((self.number_of_regions, self.number_of_regions),
                    self.travel_transmission_multiplier, dtype=float)
        np.fill_diagonal(self.travel_transmission_multiplier, 0.0)

        self.baseline_agents_travelling_matrix = vector_world.travel_matrix
        self.contacts_matrix                   = vector_world.contacts_matrix
        self.contact_ticks_matrix              = vector_world.contact_ticks_matrix

        # Rescale
        self.baseline_agents_travelling_matrix =\
            (self.scale_factor * self.baseline_agents_travelling_matrix).astype(int)

        for r in range(self.number_of_regions):
            assert self.baseline_agents_travelling_matrix[r][r] == 0
            assert self.contacts_matrix[r][r] == 0
            assert self.contact_ticks_matrix[r][r] == 0
            assert self.travel_transmission_multiplier[r][r] == 0.0

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        vector_region.border_closure_intervention = -1
        vector_region.current_region = np.full((vector_region.number_of_agents),
                                                vector_region.id, dtype=int)

    def initial_conditions(self, sim):
        """Initial regional mixing"""

        # If using ctypes, flatten arrays accordingly
        if self.enable_ctypes:
            self.baseline_agents_travelling_matrix =\
                self.baseline_agents_travelling_matrix.flatten()
            self.contacts_matrix =\
                self.contacts_matrix.flatten()
            self.contact_ticks_matrix =\
                self.contact_ticks_matrix.flatten()
            self.travel_transmission_multiplier =\
                self.travel_transmission_multiplier.flatten()

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

    def _determine_travellers(self, vector_region_batch, agents_travelling_matrix):
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

    def _calculate_average_f(self, vector_region_batch, sum_f_by_strain, average_f,
                             facemask_transmission_multiplier):
        """Calculate average transmission force for a batch of regions"""

        for vector_region in vector_region_batch:
            self.calculate_average_f(
                c_int(vector_region.number_of_agents),
                c_int(self.number_of_strains),
                c_int(vector_region.id),
                c_double(facemask_transmission_multiplier),
                c_double(vector_region.current_region_transmission_multiplier),
                c_void_p(vector_region.current_region.ctypes.data),
                c_void_p(vector_region.current_infectiousness.ctypes.data),
                c_void_p(vector_region.current_strain.ctypes.data),
                c_void_p(vector_region.current_facemask.ctypes.data),
                c_void_p(sum_f_by_strain.ctypes.data),
                c_void_p(average_f.ctypes.data),
                c_void_p(vector_region.random_state.ctypes.data),
                pointer(vector_region.random_p)
            )

    def _travel(self, vector_region_batch, transmission_force, sum_f_by_strain, mutation_matrix,
                facemask_transmission_multiplier):
        """Calculate who gets infected for a batch of regions"""

        for vector_region in vector_region_batch:
            self.travel(
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

    def dynamics_c(self, sim, day, ticks_in_day,
                   facemask_transmission_multiplier,
                   mutation_matrix, enable_parallel,
                   num_jobs, vector_region_batches):
        """Changes to regional mixing"""

        # Reset record of who is travelling abroad
        for vector_region in sim.vector_regions:
            vector_region.current_region = np.full((vector_region.number_of_agents),
                                                   vector_region.id, dtype=int)

        # Apply border closure intervention suppressing routes if necessary
        route_suppressed = np.zeros((self.number_of_regions * self.number_of_regions), dtype=int)
        for vector_region in sim.vector_regions:
            self.suppress_routes(
                c_int(self.number_of_regions),
                c_int(vector_region.id),
                c_int(vector_region.border_closure_intervention),
                c_void_p(route_suppressed.ctypes.data)
            )
        agents_travelling_matrix = copy.deepcopy(self.baseline_agents_travelling_matrix)
        self.close_borders(
            c_int(self.number_of_regions),
            c_double(self.border_closure_factor),
            c_void_p(route_suppressed.ctypes.data),
            c_void_p(agents_travelling_matrix.ctypes.data),
            c_void_p(self.baseline_agents_travelling_matrix.ctypes.data)
        )

        # Update record of who is travelling abroad
        if enable_parallel:
            Parallel(n_jobs=num_jobs, backend="threading",
                     verbose=0)(delayed(self._determine_travellers)(vector_region_batch,
                     agents_travelling_matrix) for vector_region_batch in vector_region_batches)
        else:
            self._determine_travellers(sim.vector_regions, agents_travelling_matrix)

        # Calculate average transmission force for each region
        sum_f_by_strain = np.zeros((self.number_of_regions * self.number_of_strains), dtype=float)
        average_f = np.zeros((self.number_of_regions), dtype=float)
        self._calculate_average_f(sim.vector_regions, sum_f_by_strain, average_f,
                                  facemask_transmission_multiplier)
        transmission_force = np.zeros((self.number_of_regions *
                                       self.number_of_regions), dtype=float)
        self.average_transmission_force(
            c_int(self.number_of_regions),
            c_void_p(self.contacts_matrix.ctypes.data),
            c_void_p(self.travel_transmission_multiplier.ctypes.data),
            c_void_p(self.contact_ticks_matrix.ctypes.data),
            c_void_p(transmission_force.ctypes.data),
            c_void_p(average_f.ctypes.data),
        )

        # Calculate who gets infected
        self._travel(sim.vector_regions, transmission_force, sum_f_by_strain, mutation_matrix,
                     facemask_transmission_multiplier)

        # Infect
        for vector_region in sim.vector_regions:
            time = (day * ticks_in_day) +\
                    vector_region.prng.random_randrange(ticks_in_day)
            sim.health_model.infect_c(vector_region, time)
            sim.health_model.update_c(vector_region, day * ticks_in_day)

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
