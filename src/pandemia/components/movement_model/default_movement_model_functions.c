#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdio.h>
#include <limits.h>

// gcc -fPIC -shared -o default_movement_model_functions default_movement_model_functions.c

static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

// https://prng.di.unimi.it/xoroshiro128plus.c
uint64_t next(uint64_t * s) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = s0 + s1;
	s1 ^= s0;
	s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
	s[1] = rotl(s1, 37); // c
	return result;
}

int randrange(uint64_t * s, int num) {
    // Returns an int in the range [0, num)
    return next(s) % num;
}

int random_choice(uint64_t * s, double * weights, double sum_of_weights) {
    // Returns an int in the range [0,l) where l is the length of weights
    int rnd = randrange(s, INT_MAX);
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
    int * requesting_location_update,
    const int * requested_location_update,
    int * current_location,
    int * requesting_facemask_update,
    const int * requested_facemask_update,
    int * current_facemask
)
{
    for(int n=0; n<N; n++){
        if(requesting_location_update[n] == 1){
            current_location[n] = requested_location_update[n];
            requesting_location_update[n] = 0;
        }
        if(requesting_facemask_update[n] == 1){
            current_facemask[n] = requested_facemask_update[n];
            requesting_facemask_update[n] = 0;
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
    uint64_t * random_state
)
{
    int t_now = (t + offset) % ticks_in_week;
    int t_next = (t + 1 + offset) % ticks_in_week;

    // Determine new locations, for those agents starting a new activity
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
                    index = random_choice(random_state, weights, sum_of_weights);
                    free(weights);
                } else {
                    index = randrange(random_state, num);
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
