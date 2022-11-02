"""Default test input model, specifiying a predefined set of policy interventions"""

import logging
import numpy as np

from pandemia.components.input_model import InputModel

log = logging.getLogger("default_input_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultInputModel(InputModel):
    """Default model of input, specifiying a predefined set of policy interventions. The numbers
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

        self.lockdown_input = None
        self.border_closure_input = None
        self.facemask_input = None
        self.random_testing_input = None
        self.symptomatic_testing_input = None
        self.contact_testing_input = None
        self.vaccination_input = None

    def new_input(self, policy):
        """Set new input"""

        if policy is not None:
            self._validate_input(policy)
            self.lockdown_input            = policy.lockdown_input
            self.border_closure_input      = policy.border_closure_input
            self.facemask_input            = policy.facemask_input
            self.random_testing_input      = policy.random_testing_input
            self.symptomatic_testing_input = policy.symptomatic_testing_input
            self.contact_testing_input     = policy.contact_testing_input
            self.vaccination_input         = policy.vaccination_input
        else:
            self.lockdown_input =\
                np.full((self.simulation_length_days, self.number_of_regions), -1, dtype=int)
            self.border_closure_input =\
                np.full((self.simulation_length_days, self.number_of_regions), 1.0, dtype=float)
            self.facemask_input =\
                np.full((self.simulation_length_days, self.number_of_regions), -1, dtype=int)
            self.random_testing_input =\
                np.full((self.simulation_length_days, self.number_of_regions), 100, dtype=int)
            self.symptomatic_testing_input =\
                np.full((self.simulation_length_days, self.number_of_regions), 100, dtype=int)
            self.contact_testing_input =\
                np.full((self.simulation_length_days, self.number_of_regions), 100, dtype=int)
            self.vaccination_input =\
                np.full((self.simulation_length_days, self.number_of_regions,
                         self.number_of_vaccination_age_groups,
                         self.number_of_vaccines), 100, dtype=int)

        # Rescale
        self.random_testing_input =\
            (self.scale_factor * self.random_testing_input).astype(int)
        self.symptomatic_testing_input =\
            (self.scale_factor * self.symptomatic_testing_input).astype(int)
        self.contact_testing_input =\
            (self.scale_factor * self.contact_testing_input).astype(int)
        self.vaccination_input =\
            (self.scale_factor * self.vaccination_input).astype(int)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        pass

    def initial_conditions(self, vector_region):
        """Initial input"""

        vector_region.lockdown_intervention =\
            self.lockdown_input[0][vector_region.id]
        vector_region.facemask_intervention =\
            self.facemask_input[0][vector_region.id]
        vector_region.current_border_closure_multiplier =\
            self.border_closure_input[0][vector_region.id]
        vector_region.num_to_test_random =\
            self.random_testing_input[0][vector_region.id]
        vector_region.num_to_test_symptomatic =\
            self.symptomatic_testing_input[0][vector_region.id]
        vector_region.num_to_test_contact =\
            self.contact_testing_input[0][vector_region.id]
        vector_region.num_to_vaccinate =\
            self.vaccination_input[0][vector_region.id]

    def dynamics(self, vector_region, day):
        """Changes to input"""

        vector_region.lockdown_intervention =\
            self.lockdown_input[day][vector_region.id]
        vector_region.facemask_intervention =\
            self.facemask_input[day][vector_region.id]
        vector_region.current_border_closure_multiplier =\
            self.border_closure_input[day][vector_region.id]
        vector_region.num_to_test_random =\
            self.random_testing_input[day][vector_region.id]
        vector_region.num_to_test_symptomatic =\
            self.symptomatic_testing_input[day][vector_region.id]
        vector_region.num_to_test_contact =\
            self.contact_testing_input[day][vector_region.id]
        vector_region.num_to_vaccinate =\
            self.vaccination_input[day][vector_region.id]

    def _validate_input(self, policy):
        """Validates shape of input array"""

        assert policy.lockdown_input.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.border_closure_input.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.facemask_input.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.random_testing_input.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.symptomatic_testing_input.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.contact_testing_input.shape ==\
            (self.simulation_length_days, self.number_of_regions)
        assert policy.vaccination_input.shape ==\
            (self.simulation_length_days, self.number_of_regions,
             self.number_of_vaccination_age_groups, self.number_of_vaccines)
