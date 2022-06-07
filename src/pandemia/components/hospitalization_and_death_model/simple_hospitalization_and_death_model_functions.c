#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdio.h>

// gcc -fPIC -shared -o simple_hospitalization_and_death_model_functions.dll simple_hospitalization_and_death_model_functions.c

uint64_t next(uint64_t * s, int * p) {
    const uint64_t s0 = s[*p];
    uint64_t s1 = s[*p = (*p + 1) & 15];
    s1 ^= s1 << 31; // a
    s[*p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
    return s[*p] * 0x9e3779b97f4a7c13;
}

void dynamics_hospitalization_and_death
(
    int N, // number_of_agents
    int number_of_hospitals,
    int number_of_cemeteries,
    int id,
    double hospital_threshold,
    const int * current_location,
    const int * current_region,
    const double * current_disease,
    int * requesting_location_update,
    int * requested_location_update,
    int * requesting_facemask_update,
    int * requested_facemask_update,
    int * in_hospital,
    const int * hospitals,
    int * in_cemetery,
    const int * cemeteries,
    uint64_t * random_state,
    int * random_p
)
{
    // Hospitals
    if(number_of_hospitals > 0){
        for(int n=0; n<N; n++){
            if(current_region[n] == id && requesting_location_update[n] == 1){
                if(current_disease[n] > hospital_threshold){
                    if(in_hospital[n] == 0){
                        int index = (int) (next(random_state, random_p) % number_of_hospitals);
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
                        int index = (int) (next(random_state, random_p) % number_of_cemeteries);
                        requested_location_update[n] = cemeteries[index];
                        in_cemetery[n] = 1;
                    } else {
                        requested_location_update[n] = current_location[n];
                    }
                }
            }
        }
    }
    return;
}
