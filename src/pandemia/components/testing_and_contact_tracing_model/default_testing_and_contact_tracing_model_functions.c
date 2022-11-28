#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdio.h>
#include <limits.h>

// gcc -fPIC -shared -o default_testing_and_contact_tracing_model_functions default_testing_and_contact_tracing_model_functions.c

static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

// https://prng.di.unimi.it/xoroshiro128plus.c
uint64_t next(u_int64_t * s) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = s0 + s1;
	s1 ^= s0;
	s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
	s[1] = rotl(s1, 37); // c
	return result;
}

int randrange(u_int64_t * s,int num) {
    // Returns an int in the range [0, num)
    return next(s) % num;
}

int bernoulli(u_int64_t * s, double prob_success) {
    // Returns the value 1 if success, 0 otherwise
    int rnd = randrange(s, INT_MAX);
    int result = 0;
    if(rnd < (int) (prob_success * INT_MAX)){
        result = 1;
    }
    return result;
}

void random_shuffle(u_int64_t * s, u_int64_t * array, int n) {
    int i, j, tmp;
    for (i = n - 1; i > 0; i--){
        j = randrange(s, i + 1);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
}

void random_sample(u_int64_t * s, u_int64_t * sample, int n, u_int64_t * population, int N) {
    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    int u;
    while (m < n)
    {
        u = randrange(s, N - t);
        if ( u >= n - m )
        {
            t++;
        }
        else
        {
            sample[m] = population[t];
            t++; m++;
        }
    }
}

int test
(
    u_int64_t * s,
    double infectiousness,
    double test_threshold,
    double test_false_negative
)
{
    // Tests agent based on infectiousness
    int test_result;
    if(infectiousness >= test_threshold){
        test_result = bernoulli(s, 1 - test_false_negative);
    }
    return test_result;
}

void default_testing_and_contact_tracing_dynamics
(
    int N, // number_of_agents
    int day,
    int quarantine_period_days,
    int num_to_test_random,
    int num_to_test_symptomatic,
    int num_to_test_contact,
    int max_regular_contacts_to_test,
    int id,
    double symptomatic_disease_treshold,
    double test_threshold,
    double test_false_negative,
    double prob_quarantine_with_symptoms_without_test,
    double prob_quarantine_with_contact_without_test,
    const u_int64_t * num_regular_contacts_to_test,
    const u_int64_t * regular_contacts_to_test,
    const u_int64_t * current_region,
    const double * current_infectiousness,
    const double * current_disease,
    double * yesterdays_disease,
    u_int64_t * end_of_quarantine_days,
    u_int64_t * current_quarantine,
    u_int64_t * random_state
)
{
    // End quarantine if necessary
    for(int n=0; n<N; n++){
        if(end_of_quarantine_days[n] == day){
            current_quarantine[n] = 0;
        }
    }

    // Population eligible to be tested today
    u_int64_t * eligible = (u_int64_t *)malloc(sizeof(u_int64_t) * N);
    for(int n=0; n<N; n++){
        if(current_region[n] == id && current_disease[n] < 1.0){
            eligible[n] = 1;
        } else {
            eligible[n] = 0;
        }
    }

    // A record of agents newly testing positive today
    u_int64_t * newly_testing_positive = (u_int64_t *)malloc(sizeof(u_int64_t) * N);
    for(int n=0; n<N; n++){newly_testing_positive[n] = 0;}

    // Random testing (test eligible agents at random)
    if(num_to_test_random > 0){
        int num_eligible_rand = 0;
        for(int n=0; n<N; n++){if(eligible[n] == 1){num_eligible_rand += 1;}}
        u_int64_t * eligible_agents_rand = (u_int64_t *)malloc(sizeof(u_int64_t) * num_eligible_rand);
        int j = 0;
        for(int n=0; n<N; n++){if(eligible[n] == 1){eligible_agents_rand[j] = n; j += 1;}}
        int num_agents_to_test_random;
        num_agents_to_test_random = fmin(num_to_test_random, num_eligible_rand);
        u_int64_t * agents_to_test_random = (u_int64_t *)malloc(sizeof(u_int64_t) * num_agents_to_test_random);
        random_sample(random_state, agents_to_test_random, num_agents_to_test_random,
                      eligible_agents_rand, num_eligible_rand);
        for(int j=0; j<num_agents_to_test_random; j++){
            int n, test_result;
            n = agents_to_test_random[j];
            eligible[n] = 0;
            test_result = test(random_state, current_infectiousness[n], test_threshold,
                               test_false_negative);
            if(test_result == 1){
                newly_testing_positive[n] = 1;
                current_quarantine[n] = 1;
                end_of_quarantine_days[n] = day + quarantine_period_days;
            }
        }
        free(agents_to_test_random);
        free(eligible_agents_rand);
    }

    // Symptomatic testing (test eligible agents who have just become symptomatic)
    int num_eligible_symp = 0; for(int n=0; n<N; n++){if(eligible[n] == 1){num_eligible_symp += 1;}}
    u_int64_t * eligible_agents_symp = (u_int64_t *)malloc(sizeof(u_int64_t) * num_eligible_symp);
    int j = 0; for(int n=0; n<N; n++){if(eligible[n] == 1){eligible_agents_symp[j] = n; j += 1;}}
    random_shuffle(random_state, eligible_agents_symp, num_eligible_symp);
    for(int j=0; j<num_eligible_symp; j++){
        int n, test_result;
        n = eligible_agents_symp[j];
        if(current_disease[n] >= symptomatic_disease_treshold){
            if(yesterdays_disease[n] < symptomatic_disease_treshold){
                if(num_to_test_symptomatic > 0){
                    test_result = test(random_state, current_infectiousness[n], test_threshold,
                                       test_false_negative);
                    if(test_result == 1){
                        newly_testing_positive[n] = 1;
                        current_quarantine[n] = 1;
                        end_of_quarantine_days[n] = day + quarantine_period_days;
                    }
                    num_to_test_symptomatic -= 1;
                } else {
                    if(bernoulli(random_state, prob_quarantine_with_symptoms_without_test)){
                        current_quarantine[n] = 1;
                        end_of_quarantine_days[n] = day + quarantine_period_days;
                    }
                }
            }
        }
    }
    free(eligible_agents_symp);

    // Contact tracing (contact trace and test eligible agents)
    int num_newly_testing_positive = 0;
    for(int n=0; n<N; n++){if(newly_testing_positive[n] == 1){num_newly_testing_positive += 1;}}
    if(num_newly_testing_positive > 0 && num_to_test_contact > 0){
        u_int64_t * at_risk = (u_int64_t *)malloc(sizeof(u_int64_t) * N);
        for(int n=0; n<N; n++){at_risk[n] = 0;}
        for(int n1=0; n1<N; n1++){
            if(newly_testing_positive[n1] == 1){
                for(int j=0; j<num_regular_contacts_to_test[n1]; j++){
                    int n2;
                    n2 = regular_contacts_to_test[(n1 * max_regular_contacts_to_test) + j];
                    if(current_region[n2] == id && current_disease[n2] < 1.0){
                        at_risk[n2] = 1;
                    }
                }
            }
        }
        for(int n=0; n<N; n++){if(newly_testing_positive[n] == 1){at_risk[n] = 0;}}
        int num_agents_at_risk = 0;
        for(int n=0; n<N; n++){if(at_risk[n] == 1){num_agents_at_risk += 1;}}
        u_int64_t * agents_at_risk = (u_int64_t *)malloc(sizeof(u_int64_t) * num_agents_at_risk);
        int j = 0; for(int n=0; n<N; n++){if(at_risk[n] == 1){agents_at_risk[j] = n; j += 1;}}
        random_shuffle(random_state, agents_at_risk, num_agents_at_risk);
        for(int j=0; j<fmin(num_to_test_contact, num_agents_at_risk); j++){
            int n, test_result;
            n = agents_at_risk[j];
            test_result = test(random_state, current_infectiousness[n], test_threshold,
                               test_false_negative);
            if(test_result == 1){
                current_quarantine[n] = 1;
                end_of_quarantine_days[n] = day + quarantine_period_days;
            }
        }
        for(int j=fmin(num_to_test_contact, num_agents_at_risk); j<num_agents_at_risk; j++){
            int n;
            n = agents_at_risk[j];
            if(bernoulli(random_state, prob_quarantine_with_contact_without_test)){
                current_quarantine[n] = 1;
                end_of_quarantine_days[n] = day + quarantine_period_days;
            }
        }
        free(at_risk);
        free(agents_at_risk);
    }
    free(eligible);
    free(newly_testing_positive);

    // Update yesterdays disease
    for(int n=0; n<N; n++){
        yesterdays_disease[n] = current_disease[n];
    }
    return;
}
