#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdio.h>

// gcc -fPIC -shared -o simple_testing_and_contact_tracing_model_functions.dll simple_testing_and_contact_tracing_model_functions.c

uint64_t next(uint64_t * s, int * p) {
    const uint64_t s0 = s[*p];
    uint64_t s1 = s[*p = (*p + 1) & 15];
    s1 ^= s1 << 31; // a
    s[*p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
    return s[*p] * 0x9e3779b97f4a7c13;
}

int randrange(uint64_t * s, int * p, int num) {
    // Returns an int in the range [0, num)
    return next(s, p) % num;
}

int bernoulli(uint64_t * s, int * p, double prob_success) {
    // Returns the value 1 if success, 0 otherwise
    int rnd = randrange(s, p, INT_MAX);
    int result = 0;
    if(rnd < (int) (prob_success * INT_MAX)){
        result = 1;
    }
    return result;
}

void random_shuffle(uint64_t * s, int * p, int * array, int n) {
    int i, j, tmp;
    for (i = n - 1; i > 0; i--){
        j = randrange(s, p, i + 1);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
}

void random_sample(uint64_t * s, int * p, int * sample, int n, int * population, int N) {
    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    int u;
    while (m < n)
    {
        u = randrange(s, p, N - t);
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
    uint64_t * s,
    int * p,
    double infectiousness,
    double test_threshold,
    double test_false_negative
)
{
    // Tests agent based on infectiousness
    int test_result;
    if(infectiousness >= test_threshold){
        test_result = bernoulli(s, p, 1 - test_false_negative);
    }
    return test_result;
}

void simple_testing_and_contact_tracing_dynamics
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
    double prob_self_isolate_with_symptoms_without_test,
    double prob_quarantine_with_contact_without_test,
    const int * num_regular_contacts_to_test,
    const int * regular_contacts_to_test,
    const int * current_region,
    const double * current_infectiousness,
    const double * current_disease,
    double * yesterdays_disease,
    int * end_of_quarantine_days,
    int * current_quarantine,
    uint64_t * random_state,
    int * random_p
)
{
    // End quarantine if necessary
    for(int n=0; n<N; n++){
        if(end_of_quarantine_days[n] == day){
            current_quarantine[n] = 0;
        }
    }

    // Population eligible to be tested today
    int * eligible = (int *)malloc(sizeof(int) * N);
    for(int n=0; n<N; n++){
        if(current_region[n] == id && current_disease[n] < 1.0){
            eligible[n] = 1;
        } else {
            eligible[n] = 0;
        }
    }

    // A record of agents newly testing positive today
    int * newly_testing_positive = (int *)malloc(sizeof(int) * N);
    for(int n=0; n<N; n++){newly_testing_positive[n] = 0;}

    // Random testing
    if(num_to_test_random > 0){
        int num_eligible_rand = 0;
        for(int n=0; n<N; n++){if(eligible[n] == 1){num_eligible_rand += 1;}}
        int * eligible_agents_rand = (int *)malloc(sizeof(int) * num_eligible_rand);
        int j = 0;
        for(int n=0; n<N; n++){if(eligible[n] == 1){eligible_agents_rand[j] = n; j += 1;}}
        int num_agents_to_test_random;
        num_agents_to_test_random = fmin(num_to_test_random, num_eligible_rand);
        int * agents_to_test_random = (int *)malloc(sizeof(int) * num_agents_to_test_random);
        random_sample(random_state, random_p, agents_to_test_random,
                        num_agents_to_test_random, eligible_agents_rand, num_eligible_rand);
        for(int j=0; j<num_agents_to_test_random; j++){
            int n, test_result;
            n = agents_to_test_random[j];
            eligible[n] = 0;
            test_result = test(random_state, random_p, current_infectiousness[n],
                               test_threshold, test_false_negative);
            if(test_result == 1){
                newly_testing_positive[n] = 1;
                current_quarantine[n] = 1;
                end_of_quarantine_days[n] = day + quarantine_period_days;
            }
        }
        free(agents_to_test_random);
        free(eligible_agents_rand);
    }

    // Symptomatic testing
    int num_eligible_symp = 0; for(int n=0; n<N; n++){if(eligible[n] == 1){num_eligible_symp += 1;}}
    int * eligible_agents_symp = (int *)malloc(sizeof(int) * num_eligible_symp);
    int j = 0; for(int n=0; n<N; n++){if(eligible[n] == 1){eligible_agents_symp[j] = n; j += 1;}}
    random_shuffle(random_state, random_p, eligible_agents_symp, num_eligible_symp);
    for(int j=0; j<num_eligible_symp; j++){
        int n, test_result;
        n = eligible_agents_symp[j];
        if(current_disease[n] >= symptomatic_disease_treshold){
            if(yesterdays_disease[n] < symptomatic_disease_treshold){
                if(num_to_test_symptomatic > 0){
                    test_result = test(random_state, random_p, current_infectiousness[n],
                                       test_threshold, test_false_negative);
                    if(test_result == 1){
                        newly_testing_positive[n] = 1;
                        current_quarantine[n] = 1;
                        end_of_quarantine_days[n] = day + quarantine_period_days;
                    }
                    num_to_test_symptomatic -= 1;
                } else {
                    if(bernoulli(random_state, random_p,
                                 prob_self_isolate_with_symptoms_without_test)){
                        current_quarantine[n] = 1;
                        end_of_quarantine_days[n] = day + quarantine_period_days;
                    }
                }
            }
        }
    }
    free(eligible_agents_symp);

    // Contact tracing
    int num_newly_testing_positive = 0;
    for(int n=0; n<N; n++){if(newly_testing_positive[n] == 1){num_newly_testing_positive += 1;}}
    if(num_newly_testing_positive > 0 && num_to_test_contact > 0){
        int * at_risk = (int *)malloc(sizeof(int) * N);
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
        int * agents_at_risk = (int *)malloc(sizeof(int) * num_agents_at_risk);
        int j = 0; for(int n=0; n<N; n++){if(at_risk[n] == 1){agents_at_risk[j] = n; j += 1;}}
        random_shuffle(random_state, random_p, agents_at_risk, num_agents_at_risk);
        for(int j=0; j<fmin(num_to_test_contact, num_agents_at_risk); j++){
            int n, test_result;
            n = agents_at_risk[j];
            test_result = test(random_state, random_p, current_infectiousness[n],
                                test_threshold, test_false_negative);
            if(test_result == 1){
                current_quarantine[n] = 1;
                end_of_quarantine_days[n] = day + quarantine_period_days;
            }
        }
        for(int j=fmin(num_to_test_contact, num_agents_at_risk); j<num_agents_at_risk; j++){
            int n;
            n = agents_at_risk[j];
            if(bernoulli(random_state, random_p, prob_quarantine_with_contact_without_test)){
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
