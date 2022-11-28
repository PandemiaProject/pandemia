#include <inttypes.h>

// gcc -fPIC -shared -o simulator_functions simulator_functions.c

void collect_telemetry_data
(
    int N, // number_of_agents
    int S, // number_of_strains
    int id,
    const uint64_t * current_strain,
    uint64_t * strain_counts
)
{
    for(int n=0; n<N; n++){
        if(current_strain[n] != -1){
            strain_counts[(id * S) + current_strain[n]] += 1;
        }
    }
    }

int count_dead
(
    int N, // number_of_agents
    const double * current_disease
)
{
    int deaths = 0;
    for(int n=0; n<N; n++){
        if(current_disease[n] == 1.0){
            deaths += 1;
        }
    }
    return deaths;
}