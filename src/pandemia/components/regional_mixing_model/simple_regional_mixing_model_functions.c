#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdio.h>

// gcc -fPIC -shared -o simple_regional_mixing_model_functions.dll simple_regional_mixing_model_functions.c

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

int random_choice(uint64_t * s, int * p, double * weights, double sum_of_weights) {
    // Returns an int in the range [0,l) where l is the length of weights
    int rnd = randrange(s, p, INT_MAX);
    int index = 0;
    while (rnd >= (int) ((weights[index] / sum_of_weights) * INT_MAX)){
        rnd -= (weights[index] / sum_of_weights) * INT_MAX;
        ++index;
    }
    return index;
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

void close_borders
(
    int R, // number_of_regions
    int id,
    double current_border_closure_multiplier,
    int * agents_travelling_matrix,
    const int * baseline_agents_travelling_matrix
)
{
    int r1 = id;
    for(int r2=0; r2<R; r2++){
        agents_travelling_matrix[(r1 * R) + r2] =
            (int) ((double) baseline_agents_travelling_matrix[(r1 * R) + r2] *
                            current_border_closure_multiplier);
    }
    return;
}

void determine_travellers
(
    int N, // number_of_agents
    int id,
    int R, // number_of_regions
    const double * current_disease,
    const int * current_strain,
    int * current_region,
    const int * agents_travelling_matrix,
    uint64_t * random_state,
    int * random_p
)
{
    int r1 = id;
    int total_num_to_travel = 0;

    // Determine how many agents travel from this region today
    for(int r2=0; r2<R; r2++){total_num_to_travel += agents_travelling_matrix[(r1 * R) + r2];}

    if(total_num_to_travel > 0){

        // Determine who is eligible to travel from this region today
        int * eligible = (int *)malloc(sizeof(int) * N);
        for(int n=0; n<N; n++){
            if(current_strain[n] == -1 && current_disease[n] < 1.0){
                eligible[n] = 1;
            } else {
                eligible[n] = 0;
            }
        }

        // Among those eligible determine who travels today
        for(int r2=0; r2<R; r2++){
            int num_eligible_to_travel = 0;
            for(int n=0; n<N; n++){if(eligible[n] == 1){num_eligible_to_travel += 1;}}
            int * agents_eligible_to_travel = (int *)malloc(sizeof(int) * num_eligible_to_travel);
            int j = 0;
            for(int n=0; n<N; n++){if(eligible[n] == 1){agents_eligible_to_travel[j] = n; j += 1;}}
            int num_to_travel;
            num_to_travel = agents_travelling_matrix[(r1 * R) + r2];
            if(num_to_travel > 0){
                int num_agents_to_travel;
                num_agents_to_travel = fmin(num_to_travel, num_eligible_to_travel);
                int * agents_to_travel = (int *)malloc(sizeof(int) * num_agents_to_travel);
                random_sample(random_state, random_p, agents_to_travel,
                              num_agents_to_travel, agents_eligible_to_travel,
                              num_eligible_to_travel);
                for(int j=0; j<num_agents_to_travel; j++){
                    int n;
                    n = agents_to_travel[j];
                    eligible[n] = 0;
                    current_region[n] = r2;
                }
                free(agents_to_travel);
            }
            free(agents_eligible_to_travel);
        }
        free(eligible);
    }
    return;
}

void transmission_out
(
    int N, // number_of_agents
    int S, // number_of_strains
    int id,
    double facemask_transmission_multiplier,
    double travel_multiplier,
    double current_region_transmission_multiplier,
    const int * current_region,
    const double * current_infectiousness,
    const int * current_strain,
    const int * current_facemask,
    double * sum_f_by_strain,
    double * transmission_force,
    uint64_t * random_state,
    int * random_p
)
{
    double facemask_multiplier, f;
    for(int n=0; n<N; n++){
        if(current_region[n] == id && current_strain[n] != -1){
            facemask_multiplier = 1 + current_facemask[n] * (facemask_transmission_multiplier - 1);
            f = current_region_transmission_multiplier *
                facemask_multiplier * current_infectiousness[n];
            f = (f * travel_multiplier) / N;
            sum_f_by_strain[(id * S) + current_strain[n]] += f;
            transmission_force[id] *= 1 - f;
        }
    }
    transmission_force[id] = 1 - transmission_force[id];
    return;
}

void transmission_in
(
    int R, // number_of_regions
    int S, // number_of_strains
    int N, // number_of_agents
    int r1,
    int * current_facemask,
    int * current_region,
    double facemask_transmission_multiplier,
    double * sum_f_by_strain,
    double * current_sigma_immunity_failure,
    int * infection_event,
    double * transmission_force,
    double * mutation_matrix,
    uint64_t * random_state,
    int * random_p
)
{
    double * weights = (double *)malloc(sizeof(double) * S);
    double sum_of_weights;
    int s1, s2;
    for(int n=0; n<N; n++){
        if(current_region[n] != r1){
            int r2 = current_region[n];
            double facemask_multiplier;
            facemask_multiplier =
                1 + current_facemask[n] * (facemask_transmission_multiplier - 1);
            double prob;
            if(S > 1){
                prob = facemask_multiplier * transmission_force[r2];
                if(bernoulli(random_state, random_p, prob) == 1){
                    sum_of_weights = 0;
                    for(int s=0; s<S; s++){
                        weights[s] = sum_f_by_strain[(r2 * S) + s];
                        sum_of_weights += weights[s];
                    }
                    s1 = random_choice(random_state, random_p, weights, sum_of_weights);
                    double prob;
                    prob = current_sigma_immunity_failure[(n * S) + s1];
                    if(bernoulli(random_state, random_p, prob) == 1){
                        sum_of_weights = 0;
                        for(int s=0; s<S; s++){
                            weights[s] = mutation_matrix[(s1 * S) + s];
                            sum_of_weights += weights[s];
                        }
                        s2 = random_choice(random_state, random_p, weights, sum_of_weights);
                        infection_event[n] = s2;
                    }
                }
            } else {
                s1 = 0;
                prob = facemask_multiplier *
                       transmission_force[r2] *
                       current_sigma_immunity_failure[(n * S) + s1];
                if(bernoulli(random_state, random_p, prob) == 1){
                    infection_event[n] = s1;
                }
            }
        }

    }
    free(weights);
    return;
}
