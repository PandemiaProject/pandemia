#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdio.h>

// gcc -fPIC -shared -o default_hospitalization_and_death_model_functions default_hospitalization_and_death_model_functions.c

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

void dynamics_hospitalization_and_death
(
    int N, // number_of_agents
    int number_of_hospitals,
    int number_of_cemeteries,
    int id,
    double hospital_threshold,
    const uint64_t * current_location,
    const uint64_t * current_region,
    const double * current_disease,
    const uint64_t * requesting_location_update,
    uint64_t * requested_location_update,
    uint64_t * requesting_facemask_update,
    uint64_t * requested_facemask_update,
    uint64_t * in_hospital,
    const uint64_t * hospitals,
    uint64_t * in_cemetery,
    const uint64_t * cemeteries,
    uint64_t * random_state
)
{
    // Hospitals
    if(number_of_hospitals > 0){
        for(int n=0; n<N; n++){
            if(current_region[n] == id && requesting_location_update[n] == 1){
                if(current_disease[n] > hospital_threshold){
                    if(in_hospital[n] == 0){
                        int index = (int) (next(random_state) % number_of_hospitals);
                        requested_location_update[n] = hospitals[index];
                        in_hospital[n] = 1;
                    } else {
                        requested_location_update[n] = current_location[n];
                    }
                    requesting_facemask_update[n] = 1;
                    requested_facemask_update[n] = 1;
                } else {
                    if(in_hospital[n] == 1){
                        in_hospital[n] = 0;
                    }
                }
            }
        }
    }

    // Cemeteries
    if(number_of_cemeteries > 0){
        for(int n=0; n<N; n++){
            if(current_region[n] == id && requesting_location_update[n] == 1){
                if(current_disease[n] == 1.0){
                    if(in_cemetery[n] == 0){
                        int index = (int) (next(random_state) % number_of_cemeteries);
                        requested_location_update[n] = cemeteries[index];
                        in_cemetery[n] = 1;
                    } else {
                        requested_location_update[n] = current_location[n];
                    }
                }
            }
        }
    }
    }
