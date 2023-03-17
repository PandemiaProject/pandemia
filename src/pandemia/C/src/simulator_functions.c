#include <inttypes.h>

// gcc -fPIC -shared -o simulator_functions simulator_functions.c

void collect_telemetry_data
(
    int32_t N, // number_of_agents
    int32_t S, // number_of_strains
    int32_t id,
    const int32_t * current_strain,
    int64_t * infections
)
{
    for(int32_t n=0; n<N; n++){
        if(current_strain[n] != -1){
            infections[(id * S) + current_strain[n]] += 1;
        }
    }
    }

int32_t count_dead
(
    int32_t N, // number_of_agents
    const double * current_disease
)
{
    int32_t deaths = 0;
    for(int32_t n=0; n<N; n++){
        if(current_disease[n] == 1.0){
            deaths += 1;
        }
    }
    return deaths;
}