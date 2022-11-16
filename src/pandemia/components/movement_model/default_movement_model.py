"""Default random location selection."""

import logging
import numpy as np

from ctypes import c_void_p, c_int, cdll

from ..movement_model import MovementModel

log = logging.getLogger("default_movement_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultMovementModel(MovementModel):
    """Uses random sampling to select locations in response to activity changes. When an agent
    n selects a new activity a, a new location is chosen for them at random from the array
    activity_locations[n][a]. If weighted sampling is used, then the random choice is made using
    the weights contained in the array activity_location_weights[n][a]. Location closures are also
    implemented here, which can be used to represent, for example, lockdowns. The wearing of
    facemasks is also determined here, since the wearing of facemasks is assumed to depend on the
    location type and activity of agents."""

    def __init__(self, config):
        """Initial agent locations"""
        super().__init__(config)

        lib = cdll.LoadLibrary("./src/pandemia/components/movement_model/"
                                "default_movement_model_functions")

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
        """Initialize location of each agent"""

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

        # Determine which locations may be subject to closure
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

        # Flatten arrays
        vector_region.weekly_routines = vector_region.weekly_routines.flatten()
        vector_region.wears_facemask = vector_region.wears_facemask.flatten()
        vector_region.activity_locations = vector_region.activity_locations.flatten()
        vector_region.activity_location_weights =\
            vector_region.activity_location_weights.flatten()
        vector_region.num_activity_locations = vector_region.num_activity_locations.flatten()
        vector_region.location_closure = vector_region.location_closure.flatten()

    def update(self, vector_region, t):
        """Updates related to movement"""

        self.update_movement(
            c_int(vector_region.number_of_agents),
            c_void_p(vector_region.requesting_location_update.ctypes.data),
            c_void_p(vector_region.requested_location_update.ctypes.data),
            c_void_p(vector_region.current_location.ctypes.data),
            c_void_p(vector_region.requesting_facemask_update.ctypes.data),
            c_void_p(vector_region.requested_facemask_update.ctypes.data),
            c_void_p(vector_region.current_facemask.ctypes.data)
        )

    def dynamics(self, vector_region, t, offset, ticks_in_week):
        """Changes related to movement"""

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
            c_void_p(vector_region.random_state.ctypes.data)
        )

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
