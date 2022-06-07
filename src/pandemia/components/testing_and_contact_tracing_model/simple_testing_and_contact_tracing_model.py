"""Simple testing and contact tracing model"""

from collections import defaultdict
import logging
import numpy as np

from ctypes import c_void_p, c_double, pointer, c_int, cdll

from pandemia.components.testing_and_contact_tracing_model import TestingAndContactTracingModel

log = logging.getLogger("simple_testing_and_contact_tracing_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class SimpleTestingAndContactTracingModel(TestingAndContactTracingModel):
    """Simple model of testing and contact tracing"""

    def __init__(self, config, enable_ctypes):
        """Initialize component"""
        super().__init__(config)

        self.enable_ctypes = enable_ctypes

        if self.enable_ctypes:

            lib = cdll.LoadLibrary("./src/pandemia/components/testing_and_contact_tracing_model/"
                                   "simple_testing_and_contact_tracing_model_functions.dll")

            self.simple_testing_and_contact_tracing_dynamics =\
                lib.simple_testing_and_contact_tracing_dynamics

            self.simple_testing_and_contact_tracing_dynamics.restype = None

        self.quarantine_period_days       = config['quarantine_period_days']
        self.symptomatic_disease_treshold = config['symptomatic_disease_treshold']
        self.test_threshold               = config['test_threshold']
        self.test_false_negative          = config['test_false_negative']
        self.max_regular_contacts_to_test = config['max_regular_contacts_to_test']

        self.prob_self_isolate_with_symptoms_without_test =\
            config['prob_self_isolate_with_symptoms_without_test']
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

        # If using ctypes, flatten arrays accordingly
        if self.enable_ctypes:
            vector_region.regular_contacts_to_test =\
                vector_region.regular_contacts_to_test.flatten()

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
        """Change to testing and contact tracing"""

        if self.enable_ctypes:
            self.dynamics_c(vector_region, day)
        else:
            self.dynamics_python(vector_region, day)

    def dynamics_c(self, vector_region, day):
        """Change to testing and contact tracing"""

        self.simple_testing_and_contact_tracing_dynamics(
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
            c_double(self.prob_self_isolate_with_symptoms_without_test),
            c_double(self.prob_quarantine_with_contact_without_test),
            c_void_p(vector_region.num_regular_contacts_to_test.ctypes.data),
            c_void_p(vector_region.regular_contacts_to_test.ctypes.data),
            c_void_p(vector_region.current_region.ctypes.data),
            c_void_p(vector_region.current_infectiousness.ctypes.data),
            c_void_p(vector_region.current_disease.ctypes.data),
            c_void_p(vector_region.yesterdays_disease.ctypes.data),
            c_void_p(vector_region.end_of_quarantine_days.ctypes.data),
            c_void_p(vector_region.current_quarantine.ctypes.data),
            c_void_p(vector_region.random_state.ctypes.data),
            pointer(vector_region.random_p)
        )

    def dynamics_python(self, vector_region, day):
        """Change to testing and contact tracing"""

        # End quarantine if necessary
        for n in range(vector_region.number_of_agents):
            if vector_region.end_of_quarantine_days[n] == day:
                vector_region.current_quarantine[n] = 0

        # Population eligible to be tested today
        eligible = [n for n in range(vector_region.number_of_agents)
                    if vector_region.current_region[n] == vector_region.id and
                    vector_region.current_disease[n] < 1.0]

        # A record of agents newly testing positive today
        newly_testing_positive = []

        # Random testing
        if vector_region.num_to_test_random > 0:
            agents_to_test_random = vector_region.prng.random_sample(eligible,
                                        min(vector_region.num_to_test_random, len(eligible)))
            for n in agents_to_test_random:
                eligible.remove(n)
                test_result = self._test(vector_region.current_infectiousness[n],
                                         vector_region.prng)
                if test_result == 1:
                    newly_testing_positive.append(n)
                    vector_region.current_quarantine[n] = 1
                    vector_region.end_of_quarantine_days[n] = day + self.quarantine_period_days

        # Symptomatic testing
        vector_region.prng.random_shuffle(eligible)
        for n in eligible:
            if vector_region.current_disease[n] >= self.symptomatic_disease_treshold:
                if vector_region.yesterdays_disease[n] < self.symptomatic_disease_treshold:
                    if vector_region.num_to_test_symptomatic > 0:
                        test_result = self._test(vector_region.current_infectiousness[n],
                                                 vector_region.prng)
                        if test_result == 1:
                            newly_testing_positive.append(n)
                            vector_region.current_quarantine[n] = 1
                            vector_region.end_of_quarantine_days[n] =\
                                day + self.quarantine_period_days
                        vector_region.num_to_test_symptomatic -= 1
                    else:
                        if vector_region.prng.random_float(1.0) <\
                            self.prob_self_isolate_with_symptoms_without_test:
                            vector_region.current_quarantine[n] = 1
                            vector_region.end_of_quarantine_days[n] =\
                                day + self.quarantine_period_days

        # Contact tracing
        if len(newly_testing_positive) > 0 and vector_region.num_to_test_contact > 0:
            at_risk = [0 for _ in range(vector_region.number_of_agents)]
            for n1 in newly_testing_positive:
                for n2 in [vector_region.regular_contacts_to_test[n1][j] for j
                           in range(vector_region.num_regular_contacts_to_test[n1])]:
                    if vector_region.current_region[n] == vector_region.id:
                        if vector_region.current_disease[n2] < 1.0:
                            at_risk[n2] = 1
            for n in newly_testing_positive:
                at_risk[n] = 0
            agents_at_risk = [n for n in range(vector_region.number_of_agents) if at_risk[n] == 1]
            num_agents_at_risk = len(agents_at_risk)
            vector_region.prng.random_shuffle(agents_at_risk)
            to_test = agents_at_risk[0:min(vector_region.num_to_test_contact, num_agents_at_risk)]
            not_to_test =\
                agents_at_risk[min(vector_region.num_to_test_contact, num_agents_at_risk):]
            for n in to_test:
                test_result = self._test(vector_region.current_infectiousness[n],
                                         vector_region.prng)
                if test_result == 1:
                    vector_region.current_quarantine[n] = 1
                    vector_region.end_of_quarantine_days[n] = day + self.quarantine_period_days
            for n in not_to_test:
                if vector_region.prng.random_float(1.0) <\
                    self.prob_quarantine_with_contact_without_test:
                    vector_region.current_quarantine[n] = 1
                    vector_region.end_of_quarantine_days[n] = day + self.quarantine_period_days

        # Update yesterdays disease
        for n in range(vector_region.number_of_agents):
            vector_region.yesterdays_disease[n] = vector_region.current_disease[n]

    def _test(self, infectiousness, prng):
        """Tests agent based on infectiousness"""

        test_result = 0
        if infectiousness >= self.test_threshold:
            if prng.random_float(1.0) < 1 - self.test_false_negative:
                test_result = 1

        return test_result
