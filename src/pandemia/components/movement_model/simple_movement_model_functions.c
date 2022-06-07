#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdio.h>

// gcc -fPIC -shared -o simple_movement_model_functions.dll simple_movement_model_functions.c

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

void update_movement
(
    int N, // number_of_agents
    int * l_requesting,
    const int * l_requested,
    int * l_current,
    int * f_requesting,
    const int * f_requested,
    int * f_current
)
{
    for(int n=0; n<N; n++){
        if(l_requesting[n] == 1){
            l_current[n] = l_requested[n];
            l_requesting[n] = 0;
        }
        if(f_requesting[n] == 1){
            f_current[n] = f_requested[n];
            f_requesting[n] = 0;
        }
    }
    return;
}

void dynamics_movement
(
    int N, // number_of_agents
    int A, // number_of_activities
    int lockdown_intervention,
    int facemask_intervention,
    int id,
    int use_weights,
    int t,
    int offset,
    int ticks_in_week,
    int max_num_activity_locations,
    const int * current_region,
    const uint8_t * weekly_routines,
    const int * current_facemask,
    const int * wears_facemask,
    int * requested_facemask_update,
    int * requesting_facemask_update,
    const int * activity_locations,
    const double * activity_location_weights,
    const int * num_activity_locations,
    const int * location_closure,
    int * requested_location_update,
    const int * home_location,
    const int * current_quarantine,
    int * requesting_location_update,
    uint64_t * random_state,
    int * random_p
)
{
    int t_now = (t + offset) % ticks_in_week;
    int t_next = (t + 1 + offset) % ticks_in_week;

    // Determine new locations
    for(int n=0; n<N; n++){
        if(current_region[n] == id){
            int new_activity = weekly_routines[(n * ticks_in_week) + t_next];
            if(new_activity != weekly_routines[(n * ticks_in_week) + t_now]){
                // If the agent is not currently wearing a facemask but should be, put it on
                if(current_facemask[n] == 0 && facemask_intervention == 1 &&
                   wears_facemask[(n * A) + new_activity] == 1){
                    requested_facemask_update[n] = 1;
                    requesting_facemask_update[n] = 1;
                }
                // If the agent is currently wearing a facemask but shouldn't be, take it off
                if(current_facemask[n] == 1 && (facemask_intervention == -1 ||
                   wears_facemask[(n * A) + new_activity] == 0)){
                    requested_facemask_update[n] = 0;
                    requesting_facemask_update[n] = 1;
                }
                // Randomly select an index to determine the new location
                int index;
                int num = num_activity_locations[(n * A) + new_activity];
                if(use_weights == 1){
                    double * weights = (double *)malloc(sizeof(double) * num);
                    double sum_of_weights = 0;
                    int index_0 = (n * A * max_num_activity_locations) +
                                  (new_activity * max_num_activity_locations);
                    for(int j=0; j<num; j++){
                        weights[j] = activity_location_weights[index_0 + j];
                        sum_of_weights += weights[j];
                    }
                    index = random_choice(random_state, random_p, weights, sum_of_weights);
                    free(weights);
                } else {
                    index = randrange(random_state, random_p, num);
                }
                // Override location choice if necessary
                if(current_quarantine[n] == 1 || (lockdown_intervention == 1 &&
                   location_closure[(n * A * max_num_activity_locations)
                                    + (new_activity * max_num_activity_locations)
                                    + index] == 1)){
                    requesting_location_update[n] = home_location[n];
                } else {
                    requested_location_update[n] =
                     activity_locations[(n * A * max_num_activity_locations)
                                        + (new_activity * max_num_activity_locations)
                                        + index];
                }
                requesting_location_update[n] = 1;
            }
        }
    }
    return;
}
