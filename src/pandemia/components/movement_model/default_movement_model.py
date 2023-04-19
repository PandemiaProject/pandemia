"""Default random location selection."""

import logging
import numpy as np

from ctypes import c_void_p, c_int64
from ..movement_model import MovementModel

log = logging.getLogger("default_movement_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultMovementModel(MovementModel):
    """Default model of agent movement.

    Uses random sampling to select locations in response to activity changes. When an agent
    n selects a new activity a, a new location is chosen for them at random from the array
    activity_locations[n][a]. If weighted sampling is used, then the random choice is made using
    the weights contained in the array activity_location_weights[n][a]. Location closures are also
    implemented here, which can be used to represent, for example, lockdowns. The wearing of
    facemasks is also determined here, since the wearing of facemasks is assumed to depend on the
    location type and activity of agents.

    Parameters:
    -----------
    config : Config
        A Pandemia Config object. A sub-config of the full config, containing the data used to
        configure this component.
    """

    def __init__(self, config):
        """Initialize agent movement model."""
        super().__init__(config)

        self.update_movement   = self.lib.update_movement
        self.dynamics_movement = self.lib.dynamics_movement

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
        """Initializes numpy arrays associated to this component."""

        number_of_agents = vector_region.number_of_agents
        number_of_activities = vector_region.number_of_activities
        max_num_activity_locations = vector_region.max_num_activity_locations

        vector_region.current_location = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.requested_location_update = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.requesting_location_update = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.location_closure = np.ones((number_of_agents, number_of_activities,
                                          max_num_activity_locations), dtype=np.int64)
        vector_region.lockdown_intervention = 0
        vector_region.current_facemask = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.requested_facemask_update = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.requesting_facemask_update = np.zeros((number_of_agents), dtype=np.int64)
        vector_region.wears_facemask = np.zeros((number_of_agents, number_of_activities),
                                                dtype=np.int64)
        vector_region.facemask_intervention = 0
        vector_region.current_quarantine = np.zeros((number_of_agents), dtype=np.int64)

        vector_region.home_location = np.zeros((number_of_agents), dtype=np.int64)
        assert self.home_activity in vector_region.activity_strings
        home_activity_id = vector_region.activity_strings.index(self.home_activity)
        vector_region.home_location =\
            vector_region.activity_locations[np.arange(number_of_agents), home_activity_id, 0]

    def initial_conditions(self, vector_region, offset):
        """Establish initial location of each agent."""

        # Assign initial locations
        initial_activities = vector_region.weekly_routines[:, offset]
        ids = np.arange(vector_region.number_of_agents)
        locs = vector_region.activity_locations[ids, initial_activities]
        nums = vector_region.num_activity_locations[ids, initial_activities]
        indices = vector_region.prng.prng_np.randint(nums)
        vector_region.current_location = locs[np.arange(vector_region.number_of_agents), indices]

        # Determine who wears facemasks
        self.facemask_hesitancy = np.asarray(self.facemask_hesitancy)
        age_group_indices = np.searchsorted(self.age_groups, vector_region.age, side='right') - 1
        age_group_indices = np.clip(age_group_indices, 0, len(self.age_groups) - 1)
        random_probabilities = vector_region.prng.prng_np.random(vector_region.number_of_agents)
        facemask_decision = random_probabilities < self.facemask_hesitancy[age_group_indices]
        facemask_activities_mask = np.isin(vector_region.activity_strings, self.facemask_activities)
        vector_region.wears_facemask[:, facemask_activities_mask] = 1
        vector_region.wears_facemask[facemask_decision, :] = 0

        # Determine which locations may be subject to closure
        location_typ_strings = np.array(vector_region.location_typ_strings)
        location_typs = location_typ_strings[vector_region.activity_locations]
        exemption_mask = np.zeros_like(location_typs, dtype=bool)
        for location_typ, activities in self.location_closure_exemptions.items():
            location_mask = location_typs == location_typ
            activity_mask = np.isin(vector_region.activity_strings, activities)[None, :, None]
            exemption_mask |= location_mask & activity_mask
        vector_region.location_closure[exemption_mask] = 0

        # Flatten arrays
        vector_region.weekly_routines = vector_region.weekly_routines.flatten()
        vector_region.wears_facemask = vector_region.wears_facemask.flatten()
        vector_region.activity_locations = vector_region.activity_locations.flatten()
        vector_region.activity_location_weights =\
            vector_region.activity_location_weights.flatten()
        vector_region.num_activity_locations = vector_region.num_activity_locations.flatten()
        vector_region.location_closure = vector_region.location_closure.flatten()

    def update(self, vector_region, t):
        """Updates related to movement."""

        self.update_movement(
            c_int64(vector_region.number_of_agents),
            c_void_p(vector_region.requesting_location_update.ctypes.data),
            c_void_p(vector_region.requested_location_update.ctypes.data),
            c_void_p(vector_region.current_location.ctypes.data),
            c_void_p(vector_region.requesting_facemask_update.ctypes.data),
            c_void_p(vector_region.requested_facemask_update.ctypes.data),
            c_void_p(vector_region.current_facemask.ctypes.data)
        )

    def dynamics(self, vector_region, t, offset, ticks_in_week):
        """Changes related to movement."""

        self.dynamics_movement(
            c_int64(vector_region.number_of_agents),
            c_int64(vector_region.number_of_activities),
            c_int64(vector_region.lockdown_intervention),
            c_int64(vector_region.facemask_intervention),
            c_int64(vector_region.id),
            c_int64(self.use_weights),
            c_int64(t),
            c_int64(offset),
            c_int64(ticks_in_week),
            c_int64(vector_region.max_num_activity_locations),
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
