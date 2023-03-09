"""Random policy model, specifiying a randomized set of policy interventions"""

import logging
import numpy as np

from ..policy_maker_model import PolicyMakerModel

log = logging.getLogger("random_policy_maker_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class RandomPolicyMakerModel(PolicyMakerModel):
    """Randomly generated policy interventions, for test purposes.
    """

    def __init__(self, config, scale_factor, clock, number_of_regions,
                 number_of_vaccines, age_groups):
        """Initialize component"""
        super().__init__(config, scale_factor)

        self.seed = config['seed']

        self.scale_factor = scale_factor

        self.simulation_length_days = clock.simulation_length_days
        self.number_of_regions = number_of_regions
        self.number_of_vaccines = number_of_vaccines
        self.number_of_vaccination_age_groups = len(age_groups)

        self.lockdown_policy = None
        self.border_closure_policy = None
        self.facemask_policy = None
        self.random_testing_policy = None
        self.symptomatic_testing_policy = None
        self.contact_testing_policy = None
        self.vaccination_policy = None

        np.random.seed(self.seed)

    def new_policy(self, policy):
        """Set new policy"""

        shape = (self.simulation_length_days, self.number_of_regions)
        shape_vac = (self.simulation_length_days, self.number_of_regions,
                     self.number_of_vaccination_age_groups, self.number_of_vaccines)

        self.lockdown_policy = np.random.randint(0, 2, shape, dtype=np.int32)
        self.border_closure_policy = np.random.random(shape).astype(float)
        self.facemask_policy = np.random.randint(0, 2, shape, dtype=np.int32)
        self.random_testing_policy = np.random.randint(0, 100000, shape, dtype=np.int32)
        self.symptomatic_testing_policy = np.random.randint(0, 100000, shape, dtype=np.int32)
        self.contact_testing_policy = np.random.randint(0, 100000, shape, dtype=np.int32)
        self.vaccination_policy = np.random.randint(0, 100000, shape_vac, dtype=np.int32)

        # Rescale
        self.random_testing_policy =\
            (self.scale_factor * self.random_testing_policy).astype(np.int32)
        self.symptomatic_testing_policy =\
            (self.scale_factor * self.symptomatic_testing_policy).astype(np.int32)
        self.contact_testing_policy =\
            (self.scale_factor * self.contact_testing_policy).astype(np.int32)
        self.vaccination_policy =\
            (self.scale_factor * self.vaccination_policy).astype(np.int32)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        pass

    def initial_conditions(self, vector_region):
        """Initial policy"""

        vector_region.lockdown_intervention =\
            self.lockdown_policy[0][vector_region.id]
        vector_region.facemask_intervention =\
            self.facemask_policy[0][vector_region.id]
        vector_region.current_border_closure_multiplier =\
            self.border_closure_policy[0][vector_region.id]
        vector_region.num_to_test_random =\
            self.random_testing_policy[0][vector_region.id]
        vector_region.num_to_test_symptomatic =\
            self.symptomatic_testing_policy[0][vector_region.id]
        vector_region.num_to_test_contact =\
            self.contact_testing_policy[0][vector_region.id]
        vector_region.num_to_vaccinate =\
            self.vaccination_policy[0][vector_region.id]

    def dynamics(self, vector_region, day):
        """Changes to policy"""

        vector_region.lockdown_intervention =\
            self.lockdown_policy[day][vector_region.id]
        vector_region.facemask_intervention =\
            self.facemask_policy[day][vector_region.id]
        vector_region.current_border_closure_multiplier =\
            self.border_closure_policy[day][vector_region.id]
        vector_region.num_to_test_random =\
            self.random_testing_policy[day][vector_region.id]
        vector_region.num_to_test_symptomatic =\
            self.symptomatic_testing_policy[day][vector_region.id]
        vector_region.num_to_test_contact =\
            self.contact_testing_policy[day][vector_region.id]
        vector_region.num_to_vaccinate =\
            self.vaccination_policy[day][vector_region.id]
