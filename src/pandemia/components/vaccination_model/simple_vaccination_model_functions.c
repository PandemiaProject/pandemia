#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdio.h>

// gcc -fPIC -shared -o simple_vaccination_model_functions.dll simple_vaccination_model_functions.c

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

double rho_evaluate(int * partition, double * values, int length, int r, int R, int t){

    double result;

    if(t < partition[1]){
        result = values[(0 * R) + r];
    } else {
        int j;
        j = length - 1;
        while(t < partition[j]){
            j = j - 1;
        }
        result = values[(j * R) + r];
    }
    return result;
}

double sigma_evaluate(int * partition, double * values, int length, int t){

    double result;

    if(t < partition[1]){
        result = values[0];
    } else {
        int j;
        j = length - 1;
        while(t < partition[j]){
            j = j - 1;
        }
        result = values[j];
    }
    return result;
}

void dynamics_vaccination
(
    int day,
    int ticks_in_day,
    int N, // number_of_agents
    int S, // number_of_strains
    int I, // immunity_length
    int R, // number_of_rho_immunity_outcomes
    int A, // number_of_age_groups
    int V, // number_of_vaccines
    int W, // max_vaccine_length_immunity
    int booster_waiting_time_days,
    int immunity_period_ticks,
    int id,
    const double * vaccine_rho_immunity_failure_values,
    const int * vaccine_rho_immunity_failure_partitions,
    const int * vaccine_rho_immunity_failure_lengths,
    const double * vaccine_sigma_immunity_failure_values,
    const int * vaccine_sigma_immunity_failure_partitions,
    const int * vaccine_sigma_immunity_failure_lengths,
    const int * num_to_vaccinate,
    double * rho_immunity_failure_values,
    double * sigma_immunity_failure_values,
    const int * current_region,
    const double * current_disease,
    const int * current_strain,
    int * requesting_immunity_update,
    int * most_recent_first_dose,
    int * vaccine_hesitant,
    uint64_t * random_state,
    int * random_p
)
{

    for(int age_group_index=0; age_group_index<A; age_group_index++){

        // Determine who is eligible to be vaccinated today for this age group
        int num_eligible = 0;
        int * eligible = (int *)malloc(sizeof(int) * N);
        for(int n=0; n<N; n++){
            if(current_region[n] == id &&
                current_disease[n] < 1.0 &&
                day - most_recent_first_dose[n] >= booster_waiting_time_days &&
                current_strain[n] == -1 &&
                vaccine_hesitant[n] == 0){
                eligible[n] = 1;
                num_eligible += 1;
            } else {
                eligible[n] = 0;
            }
        }
        int * eligible_agents = (int *)malloc(sizeof(int) * num_eligible);
        int j = 0;
        for(int n=0; n<N; n++){if(eligible[n] == 1){eligible_agents[j] = n; j += 1;}}

        // Vaccinate a sample of the eligible population and update their immunity functions
        int t = day * ticks_in_day;
        int num_agents_to_vaccinate = 0;
        for(int v=0; v<V; v++){
            num_agents_to_vaccinate += num_to_vaccinate[(age_group_index * A) + v];
        }
        num_agents_to_vaccinate = fmin(num_eligible, num_agents_to_vaccinate);
        int * agents_to_vaccinate = (int *)malloc(sizeof(int) * num_agents_to_vaccinate);
        random_sample(random_state, random_p, agents_to_vaccinate,
                      num_agents_to_vaccinate, eligible_agents, num_eligible);
        int v = 0;
        int num_vaccinated = 0;
        for(int j=0; j<num_agents_to_vaccinate; j++){
            int n;
            n = agents_to_vaccinate[j];
            if(v < V){

                // Determine new rho immunity

                int * rho_part = (int *)malloc(sizeof(int) * W);
                double * rho_values = (double *)malloc(sizeof(double) * W * R);
                int rho_length;

                for(int s=0; s<S; s++){

                    for(int i=0; i<W; i++){

                        for(int r=0; r<R; r++){

                            rho_values[(i * R) + r] =
                                vaccine_rho_immunity_failure_values[(v * S * W * R) +
                                                                    (s * W * R) +
                                                                    (i * R) +
                                                                    r];


                        }

                        rho_part[i] =
                            vaccine_rho_immunity_failure_partitions[(v * S * W) +
                                                                    (s * W) +
                                                                    i] + t;
                        rho_part[0] = -1;

                    }

                    rho_length =
                        vaccine_rho_immunity_failure_lengths[(v * S) + s];

                    for(int i=0; i<I; i++){
                        for(int r=0; r<R; r++){

                            rho_immunity_failure_values[(n * S * I * R) +
                                                        (s * I * R) +
                                                        (i * R) +
                                                        r] =
                                rho_immunity_failure_values[(n * S * I * R) +
                                                            (s * I * R) +
                                                            (i * R) +
                                                            r] *
                                rho_evaluate(rho_part, rho_values, rho_length, r, R,
                                             (i + 1) * immunity_period_ticks);

                        }
                    }
                }
                free(rho_part);
                free(rho_values);

                // Determine new sigma immunity

                int * sigma_part = (int *)malloc(sizeof(int) * W);
                double * sigma_values = (double *)malloc(sizeof(double) * W);
                int sigma_length;

                for(int s=0; s<S; s++){

                    for(int i=0; i<W; i++){

                        sigma_values[i] =
                            vaccine_sigma_immunity_failure_values[(v * S * W) +
                                                                  (s * W) +
                                                                  i];

                        sigma_part[i] =
                            vaccine_sigma_immunity_failure_partitions[(v * S * W) +
                                                                      (s * W) +
                                                                      i] + t;
                        sigma_part[0] = -1;

                    }

                    sigma_length =
                        vaccine_rho_immunity_failure_lengths[(v * S) + s];

                    for(int i=0; i<I; i++){

                        sigma_immunity_failure_values[(n * S * I) +
                                                      (s * I) +
                                                      i] =
                            sigma_immunity_failure_values[(n * S * I) +
                                                          (s * I) +
                                                          i] *
                            sigma_evaluate(&sigma_part[0], &sigma_values[0], sigma_length,
                                           (i + 1) * immunity_period_ticks);

                    }
                }
                free(sigma_part);
                free(sigma_values);

                most_recent_first_dose[n] = day;
                num_vaccinated += 1;
                requesting_immunity_update[n] = 1;

            }
            if(num_vaccinated >= num_to_vaccinate[(age_group_index * A) + v]){
                v += 1;
                num_vaccinated = 0;
            }
        }
        free(agents_to_vaccinate);
        free(eligible_agents);
        free(eligible);
    }
    return;
}