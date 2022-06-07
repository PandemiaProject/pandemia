"""Simple random location selection."""

import logging
import numpy as np

from ctypes import c_void_p, pointer, c_int, cdll

from pandemia.components.movement_model import MovementModel

log = logging.getLogger("simple_movement_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class SimpleMovementModel(MovementModel):
    """Uses simple random sampling to select locations in response to activity changes."""

    def __init__(self, config, enable_ctypes):
        """Initial agent locations"""
        super().__init__(config)

        self.enable_ctypes = enable_ctypes

        if self.enable_ctypes:

            lib = cdll.LoadLibrary("./src/pandemia/components/movement_model/"
                                   "simple_movement_model_functions.dll")

            self.update_movement   = lib.update_movement
            self.dynamics_movement = lib.dynamics_movement

            self.update_movement.restype   = None
            self.dynamics_movement.restype = None

        self.home_activity                      = config['home_activity']
        self.location_closure_exemptions        = config['location_closure_exemptions']
        self.facemask_activities                = config['facemask_activities']
        self.age_groups                         = config['age_groups']
        self.facemask_hesitancy                 = config['facemask_hesitancy']
        self.weighted_sampling                  = config['weighted_sampling']

        if self.weighted_sampling:
            self.use_weights = 1
        else:
            self.use_weights = 0

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        number_of_agents = vector_region.number_of_agents
        number_of_activities = vector_region.number_of_activities
        max_num_activity_locations = vector_region.max_num_activity_locations

        vector_region.current_location = np.zeros((number_of_agents), dtype=int)
        vector_region.requested_location_update = np.zeros((number_of_agents), dtype=int)
        vector_region.requesting_location_update = np.zeros((number_of_agents), dtype=int)
        vector_region.location_closure = np.ones((number_of_agents, number_of_activities,
                                          max_num_activity_locations), dtype=int)
        vector_region.lockdown_intervention = -1
        vector_region.current_facemask = np.zeros((number_of_agents), dtype=int)
        vector_region.requested_facemask_update = np.zeros((number_of_agents), dtype=int)
        vector_region.requesting_facemask_update = np.zeros((number_of_agents), dtype=int)
        vector_region.wears_facemask = np.zeros((number_of_agents, number_of_activities),
                                                dtype=int)
        vector_region.facemask_intervention = -1
        vector_region.current_quarantine = np.zeros((number_of_agents), dtype=int)

        vector_region.home_location = np.zeros((number_of_agents), dtype=int)
        assert self.home_activity in vector_region.activity_strings
        home_activity_id = vector_region.activity_strings.index(self.home_activity)
        for n in range(number_of_agents):
            vector_region.home_location[n] =\
                vector_region.activity_locations[n][home_activity_id][0]

    def initial_conditions(self, vector_region, offset):

        # Assign initial locations
        for n in range(vector_region.number_of_agents):
            initial_activity = vector_region.weekly_routines[n][offset]
            locs = vector_region.activity_locations[n][initial_activity]
            num = vector_region.num_activity_locations[n][initial_activity]
            index = vector_region.prng.random_randrange(num)
            vector_region.current_location[n] = locs[index]
            vector_region.current_quarantine[n] = 0

        # Determine who wears facemasks
        for n in range(vector_region.number_of_agents):
            age = vector_region.age[n]
            age_group_index = self._determine_age_group_index(age, self.age_groups)
            if vector_region.prng.random_float(1.0) < self.facemask_hesitancy[age_group_index]:
                for a in range(vector_region.number_of_activities):
                    vector_region.wears_facemask[n][a] = 0
            else:
                for a in range(vector_region.number_of_activities):
                    if vector_region.activity_strings[a] in self.facemask_activities:
                        vector_region.wears_facemask[n][a] = 1

        # Determine location closure
        for n in range(vector_region.number_of_agents):
            for a in range(vector_region.number_of_activities):
                locs = vector_region.activity_locations[n][a]
                num = vector_region.num_activity_locations[n][a]
                for index in range(num):
                    activity = vector_region.activity_strings[a]
                    location_typ = vector_region.location_typ_strings[locs[index]]
                    if location_typ in self.location_closure_exemptions:
                        if activity in self.location_closure_exemptions[location_typ]:
                            vector_region.location_closure[n][a][index] = 0

        # If using ctypes, flatten arrays accordingly
        if self.enable_ctypes:
            vector_region.weekly_routines = vector_region.weekly_routines.flatten()
            vector_region.wears_facemask = vector_region.wears_facemask.flatten()
            vector_region.activity_locations = vector_region.activity_locations.flatten()
            vector_region.activity_location_weights =\
                vector_region.activity_location_weights.flatten()
            vector_region.num_activity_locations = vector_region.num_activity_locations.flatten()
            vector_region.location_closure = vector_region.location_closure.flatten()

    def update(self, vector_region, t):
        """Updates movement functions"""

        if self.enable_ctypes:
            self.update_c(vector_region, t)
        else:
            self.update_python(vector_region, t)

    def update_c(self, vector_region, t):
        """Updates movement functions using precompiled C functions"""

        self.update_movement(
            c_int(vector_region.number_of_agents),
            c_void_p(vector_region.requesting_location_update.ctypes.data),
            c_void_p(vector_region.requested_location_update.ctypes.data),
            c_void_p(vector_region.current_location.ctypes.data),
            c_void_p(vector_region.requesting_facemask_update.ctypes.data),
            c_void_p(vector_region.requested_facemask_update.ctypes.data),
            c_void_p(vector_region.current_facemask.ctypes.data)
        )

    def update_python(self, vector_region, t):
        """Updates movement functions using only pure Python"""

        for n in range(vector_region.number_of_agents):
            # Update current location
            if vector_region.requesting_location_update[n] == 1:
                vector_region.current_location[n] = vector_region.requested_location_update[n]
                vector_region.requesting_location_update[n] = 0
            # Update facemasks
            if vector_region.requesting_facemask_update[n] == 1:
                vector_region.current_facemask[n] = vector_region.requested_facemask_update[n]
                vector_region.requesting_facemask_update[n] = 0

    def dynamics(self, vector_region, t, offset, ticks_in_week):
        """Changes agent locations"""

        if self.enable_ctypes:
            self.dynamics_c(vector_region, t, offset, ticks_in_week)
        else:
            self.dynamics_python(vector_region, t, offset, ticks_in_week)

    def dynamics_c(self, vector_region, t, offset, ticks_in_week):
        """Changes agent locations using precompiled C functions"""

        self.dynamics_movement(
            c_int(vector_region.number_of_agents),
            c_int(vector_region.number_of_activities),
            c_int(vector_region.lockdown_intervention),
            c_int(vector_region.facemask_intervention),
            c_int(vector_region.id),
            c_int(self.use_weights),
            c_int(t),
            c_int(offset),
            c_int(ticks_in_week),
            c_int(vector_region.max_num_activity_locations),
            c_void_p(vector_region.current_region.ctypes.data),
            c_void_p(vector_region.weekly_routines.ctypes.data),
            c_void_p(vector_region.current_facemask.ctypes.data),
            c_void_p(vector_region.wears_facemask.ctypes.data),
            c_void_p(vector_region.requested_facemask_update.ctypes.data),
            c_void_p(vector_region.requesting_facemask_update.ctypes.data),
            c_void_p(vector_region.activity_locations.ctypes.data),
            c_void_p(vector_region.activity_location_weights.ctypes.data),
            c_void_p(vector_region.num_activity_locations.ctypes.data),
            c_void_p(vector_region.location_closure.ctypes.data),
            c_void_p(vector_region.requested_location_update.ctypes.data),
            c_void_p(vector_region.home_location.ctypes.data),
            c_void_p(vector_region.current_quarantine.ctypes.data),
            c_void_p(vector_region.requesting_location_update.ctypes.data),
            c_void_p(vector_region.random_state.ctypes.data),
            pointer(vector_region.random_p)
        )

    def dynamics_python(self, vector_region, t, offset, ticks_in_week):
        """Changes agent locations using only pure Python"""

        t_now = (t + offset) % ticks_in_week
        t_next = (t + 1 + offset) % ticks_in_week

        # Determine new locations
        for n in range(vector_region.number_of_agents):
            if vector_region.current_region[n] == vector_region.id:
                new_activity = vector_region.weekly_routines[n][t_next]
                if new_activity != vector_region.weekly_routines[n][t_now]:
                    # If the agent is not currently wearing a facemask but should be, put it on
                    if vector_region.current_facemask[n] == 0 and\
                        vector_region.facemask_intervention == 1 and\
                        vector_region.wears_facemask[n][new_activity] == 1:
                        vector_region.requested_facemask_update[n] = 1
                        vector_region.requesting_facemask_update[n] = 1
                    # If the agent is currently wearing a facemask but shouldn't be, take it off
                    if vector_region.current_facemask[n] == 1 and\
                        (vector_region.facemask_intervention == -1 or\
                         vector_region.wears_facemask[n][new_activity] == 0):
                        vector_region.requested_facemask_update[n] = 0
                        vector_region.requesting_facemask_update[n] = 1
                    locs = vector_region.activity_locations[n][new_activity]
                    num = vector_region.num_activity_locations[n][new_activity]
                    weights = vector_region.activity_location_weights[n][new_activity][0:num]
                    indexes = list(range(num))
                    index = vector_region.prng.random_choices(indexes, weights[0:num], 1)[0]
                    # If the selected location is subject to location closure then send home
                    if vector_region.lockdown_intervention == 1 and\
                        vector_region.location_closure[n][new_activity][index] == 1:
                        vector_region.requested_location_update[n] = vector_region.home_location[n]
                    elif vector_region.current_quarantine[n] == 1:
                        vector_region.requested_location_update[n] = vector_region.home_location[n]
                    else:
                        vector_region.requested_location_update[n] = locs[index]
                    vector_region.requesting_location_update[n] = 1

    def _determine_age_group_index(self, age, age_groups):
        """Determine to which age group the agent belongs"""

        if age >= age_groups[-1]:
            return len(age_groups) - 1
        elif age < age_groups[0]:
            return 0
        else:
            i = 0
            while age < age_groups[i] or age >= age_groups[i+1]:
                i = i + 1
            return i
