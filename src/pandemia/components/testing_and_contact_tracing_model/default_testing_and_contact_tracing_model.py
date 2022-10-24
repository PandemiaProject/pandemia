"""Default testing and contact tracing model"""

from collections import defaultdict
import logging
import numpy as np

from ctypes import c_void_p, c_double, c_int, cdll

from pandemia.components.testing_and_contact_tracing_model import TestingAndContactTracingModel

log = logging.getLogger("default_testing_and_contact_tracing_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultTestingAndContactTracingModel(TestingAndContactTracingModel):
    """Default model of testing and contact tracing. Agents can be tested in this model for one of
    three reasons. First, they might be tested as a result of random testing of the population.
    Second, they might be tested after they develop symptoms. Third, they might be tested via the
    contract tracing system. Limited numbers of tests are available for each purpose. Testing and
    contact tracing indicate who is or who might be infected. Those who test positive are directed
    by this component to quarantine. For the duration of quarantine, these agents remain at home."""

    def __init__(self, config):
        """Initialize component"""
        super().__init__(config)

        lib = cdll.LoadLibrary("./src/pandemia/components/testing_and_contact_tracing_model/"
                                "default_testing_and_contact_tracing_model_functions.dll")

        self.default_testing_and_contact_tracing_dynamics =\
            lib.default_testing_and_contact_tracing_dynamics

        self.default_testing_and_contact_tracing_dynamics.restype = None

        self.quarantine_period_days       = config['quarantine_period_days']
        self.symptomatic_disease_treshold = config['symptomatic_disease_treshold']
        self.test_threshold               = config['test_threshold']
        self.test_false_negative          = config['test_false_negative']
        self.max_regular_contacts_to_test = config['max_regular_contacts_to_test']

        self.prob_quarantine_with_symptoms_without_test =\
            config['prob_quarantine_with_symptoms_without_test']
        self.prob_quarantine_with_contact_without_test =\
            config['prob_quarantine_with_contact_without_test']

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        number_of_agents = vector_region.number_of_agents
        max_regular_contacts_to_test = self.max_regular_contacts_to_test

        vector_region.num_to_test_random = 0
        vector_region.num_to_test_symptomatic = 0
        vector_region.num_to_test_contact = 0
        vector_region.end_of_quarantine_days = np.full((number_of_agents), -1, dtype=int)
        vector_region.max_regular_contacts_to_test = max_regular_contacts_to_test
        vector_region.num_regular_contacts_to_test = np.zeros((number_of_agents), dtype=int)
        vector_region.regular_contacts_to_test = np.zeros((number_of_agents,
                                                  max_regular_contacts_to_test), dtype=int)
        vector_region.yesterdays_disease = np.zeros((number_of_agents), dtype=float)

        # Flatten arrays
        vector_region.regular_contacts_to_test = vector_region.regular_contacts_to_test.flatten()

    def initial_conditions(self, vector_region):
        """Initial testing and contact tracing conditions"""

         # Initialize record of when agents should end quarantine
        for n in range(vector_region.number_of_agents):
            vector_region.end_of_quarantine_days[n] = -1

        # Determine regular contacts and subsample them to determine which of them will be eligible
        # for tracing and testing during the simulation
        residents = defaultdict(list)
        for n in range(vector_region.number_of_agents):
            residents[vector_region.home_location[n]].append(n)
        regular_contacts = [[] for _ in range(vector_region.number_of_agents)]
        for n in range(vector_region.number_of_agents):
            regular_contacts[n] = residents[vector_region.home_location[n]]
            regular_contacts[n].remove(n)
            vector_region.num_regular_contacts_to_test[n] =\
                min(len(regular_contacts[n]), vector_region.max_regular_contacts_to_test)
        regular_contacts_sample =\
            vector_region.prng.random_sample(regular_contacts,
                                             vector_region.num_regular_contacts_to_test[n])
        # Record the sample of regular contacts
        for index in range(0, vector_region.num_regular_contacts_to_test[n]):
            vector_region.regular_contacts_to_test[index] = regular_contacts_sample[index]
        # Pad the remainder of the array if necessary
        for index in range(vector_region.num_regular_contacts_to_test[n],
                           vector_region.max_regular_contacts_to_test):
            vector_region.regular_contacts_to_test[index] = -1

    def dynamics(self, vector_region, day):
        """Changes to testing and contact tracing"""

        self.default_testing_and_contact_tracing_dynamics(
            c_int(vector_region.number_of_agents),
            c_int(day),
            c_int(self.quarantine_period_days),
            c_int(vector_region.num_to_test_random),
            c_int(vector_region.num_to_test_symptomatic),
            c_int(vector_region.num_to_test_contact),
            c_int(self.max_regular_contacts_to_test),
            c_int(vector_region.id),
            c_double(self.symptomatic_disease_treshold),
            c_double(self.test_threshold),
            c_double(self.test_false_negative),
            c_double(self.prob_quarantine_with_symptoms_without_test),
            c_double(self.prob_quarantine_with_contact_without_test),
            c_void_p(vector_region.num_regular_contacts_to_test.ctypes.data),
            c_void_p(vector_region.regular_contacts_to_test.ctypes.data),
            c_void_p(vector_region.current_region.ctypes.data),
            c_void_p(vector_region.current_infectiousness.ctypes.data),
            c_void_p(vector_region.current_disease.ctypes.data),
            c_void_p(vector_region.yesterdays_disease.ctypes.data),
            c_void_p(vector_region.end_of_quarantine_days.ctypes.data),
            c_void_p(vector_region.current_quarantine.ctypes.data),
            c_void_p(vector_region.random_state.ctypes.data)
        )
