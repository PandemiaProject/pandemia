"""Random input model, specifiying a randomized set of policy interventions"""

import logging
import numpy as np

from ..input_model import InputModel

log = logging.getLogger("random_input_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class RandomInputModel(InputModel):
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

        self.lockdown_input = None
        self.border_closure_input = None
        self.facemask_input = None
        self.random_testing_input = None
        self.symptomatic_testing_input = None
        self.contact_testing_input = None
        self.vaccination_input = None

        np.random.seed(self.seed)

    def new_input(self, policy):
        """Set new input"""

        shape = (self.simulation_length_days, self.number_of_regions)
        shape_vac = (self.simulation_length_days, self.number_of_regions,
                     self.number_of_vaccination_age_groups, self.number_of_vaccines)

        self.lockdown_input = np.random.randint(-1, 1, shape, dtype=np.int64)
        self.border_closure_input = np.random.random(shape).astype(np.float64)
        self.facemask_input = np.random.randint(-1, 1, shape, dtype=np.int64)
        self.random_testing_input = np.random.randint(0, 100000, shape, dtype=np.int64)
        self.symptomatic_testing_input = np.random.randint(0, 100000, shape, dtype=np.int64)
        self.contact_testing_input = np.random.randint(0, 100000, shape, dtype=np.int64)
        self.vaccination_input = np.random.randint(0, 100000, shape_vac, dtype=np.int64)

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
