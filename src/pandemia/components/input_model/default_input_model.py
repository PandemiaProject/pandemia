"""Default test input model, specifiying a predefined set of policy interventions"""

import logging
import csv
import numpy as np
from os.path import exists

from ..input_model import InputModel

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

        self.input_data_filepath = config['input_data_filepath']

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

        self.lockdown_input =\
            np.full((self.simulation_length_days, self.number_of_regions), 0, dtype=np.int64)
        self.border_closure_input =\
            np.full((self.simulation_length_days, self.number_of_regions), 1.0, dtype=np.float64)
        self.facemask_input =\
            np.full((self.simulation_length_days, self.number_of_regions), 0, dtype=np.int64)
        self.random_testing_input =\
            np.full((self.simulation_length_days, self.number_of_regions), 0, dtype=np.int64)
        self.symptomatic_testing_input =\
            np.full((self.simulation_length_days, self.number_of_regions), 0, dtype=np.int64)
        self.contact_testing_input =\
            np.full((self.simulation_length_days, self.number_of_regions), 0, dtype=np.int64)
        self.vaccination_input =\
            np.full((self.simulation_length_days, self.number_of_regions,
                        self.number_of_vaccination_age_groups,
                        self.number_of_vaccines), 0, dtype=np.int64)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        iso2 = vector_region.name
        id = vector_region.id
        filepath = self.input_data_filepath + iso2 + ".csv"
        if exists(filepath):
            array = np.genfromtxt(filepath, skip_header = 1, delimiter=',', dtype=np.float64)
            array = array.T
            assert array.shape[1] == self.simulation_length_days
            self.lockdown_input[:,id]            = array[1].astype(np.int64)
            self.border_closure_input[:,id]      = (1 - array[2]).astype(np.float64)
            self.facemask_input[:,id]            = array[3].astype(np.int64)
            self.random_testing_input[:,id]      = (self.scale_factor * array[4]).astype(np.int64)
            self.symptomatic_testing_input[:,id] = (self.scale_factor * array[5]).astype(np.int64)
            self.contact_testing_input[:,id]     = (self.scale_factor * array[6]).astype(np.int64)
            if self.number_of_vaccines == 1:
                self.vaccination_input[:,id,0,0] = (self.scale_factor * array[7]).astype(np.int64)
                self.vaccination_input[:,id,1,0] = (self.scale_factor * array[8]).astype(np.int64)
                self.vaccination_input[:,id,2,0] = (self.scale_factor * array[9]).astype(np.int64)
            if self.number_of_vaccines == 2:
                self.vaccination_input[:,id,0,0] = (self.scale_factor * array[7]).astype(np.int64)
                self.vaccination_input[:,id,1,0] = (self.scale_factor * array[8]).astype(np.int64)
                self.vaccination_input[:,id,2,0] = (self.scale_factor * array[9]).astype(np.int64)
                self.vaccination_input[:,id,0,1] = (self.scale_factor * array[10]).astype(np.int64)
                self.vaccination_input[:,id,1,1] = (self.scale_factor * array[11]).astype(np.int64)
                self.vaccination_input[:,id,2,1] = (self.scale_factor * array[12]).astype(np.int64)

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
