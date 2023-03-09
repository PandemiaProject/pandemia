"""Default policy maker model, specifiying a predefined set of policy interventions"""

import logging
import csv
import numpy as np
from os.path import exists

from ..policy_maker_model import PolicyMakerModel

log = logging.getLogger("default_policy_maker_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultPolicyMakerModel(PolicyMakerModel):
    """Default model of a policy maker.

    Default model of a policy maker, specifiying a predefined set of policy interventions. This
    policy is specfied using arrays of integers and float, loaded from a csv file that can be
    editted. The filepath to this csv file should be found in the config.

    Parameters:
    -----------
    config : Config
        A Pandemia Config object. A sub-config of the full config, containing the data used to
        configure this component.
    scale_factor : float
        The scale factor, coming from the full config.
    clock : Clock
        A Pandemia Clock object, discretizing the day.
    number_of_regions : int
        The number of regions appearing the model.
    number_of_vaccines : int
        The number of vaccines appearing in the model.
    age_groups : list[int]
        A list of integers coming from the vaccination model, indicating the age groups for
        vaccination. For example, the list [0, 18, 65] indicates three age groups.
    """

    def __init__(self, config, scale_factor, clock, number_of_regions,
                 number_of_vaccines, age_groups):
        """Initialize component."""
        super().__init__(config, scale_factor)

        self.scale_factor = scale_factor

        self.policy_data_filepath = config['policy_data_filepath']

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
        """Sets new policy."""

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

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component."""

        iso2 = vector_region.name
        id = vector_region.id
        filepath = self.policy_data_filepath + iso2 + ".csv"
        if exists(filepath):
            array = np.genfromtxt(filepath, skip_header = 1, delimiter=',', dtype=float)
            array = array.T
            assert array.shape[1] == self.simulation_length_days
            self.lockdown_policy[:,id]            = array[1].astype(np.int32)
            self.border_closure_policy[:,id]      = (1 - array[2]).astype(float)
            self.facemask_policy[:,id]            = array[3].astype(np.int32)
            self.random_testing_policy[:,id]      = (self.scale_factor * array[4]).astype(np.int32)
            self.symptomatic_testing_policy[:,id] = (self.scale_factor * array[5]).astype(np.int32)
            self.contact_testing_policy[:,id]     = (self.scale_factor * array[6]).astype(np.int32)
            if self.number_of_vaccines == 1:
                self.vaccination_policy[:,id,0,0] = (self.scale_factor * array[7]).astype(np.int32)
                self.vaccination_policy[:,id,1,0] = (self.scale_factor * array[8]).astype(np.int32)
                self.vaccination_policy[:,id,2,0] = (self.scale_factor * array[9]).astype(np.int32)
            if self.number_of_vaccines == 2:
                self.vaccination_policy[:,id,0,0] = (self.scale_factor * array[7]).astype(np.int32)
                self.vaccination_policy[:,id,1,0] = (self.scale_factor * array[8]).astype(np.int32)
                self.vaccination_policy[:,id,2,0] = (self.scale_factor * array[9]).astype(np.int32)
                self.vaccination_policy[:,id,0,1] = (self.scale_factor * array[10]).astype(np.int32)
                self.vaccination_policy[:,id,1,1] = (self.scale_factor * array[11]).astype(np.int32)
                self.vaccination_policy[:,id,2,1] = (self.scale_factor * array[12]).astype(np.int32)

    def initial_conditions(self, vector_region):
        """Initial policy."""

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
        """Changes to policy."""

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
        """Validates shape of policy array."""

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
