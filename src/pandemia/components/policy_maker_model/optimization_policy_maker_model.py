"""Default policy maker model, specifiying a predefined set of policy interventions"""

import logging
import numpy as np

from ..policy_maker_model import PolicyMakerModel

log = logging.getLogger("default_policy_maker_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class OptimizationPolicyMakerModel(PolicyMakerModel):
    """Default model of a policy maker, specifiying a predefined set of policy interventions. The numbers
    appearing in the arrays below are for test purposes. To represents a new policy intervention,
    for example a historical set of interventions, the numbers in these arrays should be changed
    accordingly."""

    def __init__(self, config, scale_factor, clock, number_of_regions,
                 number_of_vaccines, age_groups):
        """Initialize component"""
        super().__init__(config, scale_factor)

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

    def new_policy(self, policy):
        """Set new policy"""

        if policy is not None:
            self._validate_policy(policy)
            self.lockdown_policy            = policy.lockdown_policy
            self.border_closure_policy      = policy.border_closure_policy
            self.facemask_policy            = policy.facemask_policy
            self.random_testing_policy      = policy.random_testing_policy
            self.symptomatic_testing_policy = policy.symptomatic_testing_policy
            self.contact_testing_policy     = policy.contact_testing_policy
            self.vaccination_policy         = policy.vaccination_policy
        else:
            self.lockdown_policy =\
                np.full((self.simulation_length_days, self.number_of_regions), 0, dtype=np.int32)
            self.border_closure_policy =\
                np.full((self.simulation_length_days, self.number_of_regions), 1.0, dtype=float)
            self.facemask_policy =\
                np.full((self.simulation_length_days, self.number_of_regions), 0, dtype=np.int32)
            self.random_testing_policy =\
                np.full((self.simulation_length_days, self.number_of_regions), 0, dtype=np.int32)
            self.symptomatic_testing_policy =\
                np.full((self.simulation_length_days, self.number_of_regions), 0, dtype=np.int32)
            self.contact_testing_policy =\
                np.full((self.simulation_length_days, self.number_of_regions), 0, dtype=np.int32)
            self.vaccination_policy =\
                np.full((self.simulation_length_days, self.number_of_regions,
                         self.number_of_vaccination_age_groups,
                         self.number_of_vaccines), 0, dtype=np.int32)

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

    def _validate_policy(self, policy):
        """Validates shape of policy array"""

        assert policy.lockdown_policy.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.border_closure_policy.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.facemask_policy.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.random_testing_policy.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.symptomatic_testing_policy.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.contact_testing_policy.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.vaccination_policy.shape ==\
            (self.simulation_length_days, self.number_of_regions,
             self.number_of_vaccination_age_groups, self.number_of_vaccines)
