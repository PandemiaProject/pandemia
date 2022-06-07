#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdio.h>

// gcc -fPIC -shared -o simple_health_model_functions.dll simple_health_model_functions.c

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

double rho_evaluate(int * partition, double * values, int length, int r, int R, int t){
    // Evaluates a rho immunity preset on rho immunity outcome r at time t, where R denotes the
    // total number of rho immunity outcomes
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
    // Evaluates a sigma immunity preset at time t
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

void update_health
(
    int t,
    int N, // number_of_agents
    int S, // number_of_strains
    int I, // immunity_length
    int R, // number_of_rho_immunity_outcomes
    int H, // max_preset_length_health
    int immunity_period_ticks,
    const int * i_partitions,
    int * i_indexes,
    const double * i_values,
    const int * i_lengths,
    double * i_current,
    const int * d_partitions,
    int * d_indexes,
    const double * d_values,
    const int * d_lengths,
    double * d_current,
    const int * s_partitions,
    int * s_indexes,
    const int * s_values,
    const int * s_lengths,
    int * s_current,
    const double * rho_immunity_failure_values,
    double * current_rho_immunity_failure,
    const double * sigma_immunity_failure_values,
    double * current_sigma_immunity_failure,
    int * requesting_immunity_update
)
{
    if(t % immunity_period_ticks == 0){for(int n=0; n<N; n++){requesting_immunity_update[n] = 1;}}
    int i = t / immunity_period_ticks;
    for(int n=0; n<N; n++){
        // Update infectiousness
        if(t == i_partitions[(n * H) + i_indexes[n]]){
            i_current[n] = i_values[(n * H) + i_indexes[n]];
            if (i_indexes[n] + 1 < i_lengths[n]){
                i_indexes[n] += 1;
            }
        }
        // Update disease
        if(t == d_partitions[(n * H) + d_indexes[n]]){
            d_current[n] = d_values[(n * H) + d_indexes[n]];
            if (d_indexes[n] + 1 < d_lengths[n]){
                d_indexes[n] += 1;
            }
        }
        // Update strain
        if(t == s_partitions[(n * H) + s_indexes[n]]){
            s_current[n] = s_values[(n * H) + s_indexes[n]];
            if (s_indexes[n] + 1 < s_lengths[n]){
                s_indexes[n] += 1;
            }
        }
        // Update rho and sigma immunity
        if(requesting_immunity_update[n] == 1){
            int nSR, nSIR, nS, nSI, iR;
            iR = i * R;
            nSR = n * S * R;
            nSIR = n * S * I * R;
            nS = n * S;
            nSI = n * S * I;
            for(int s=0; s<S; s++){
                int sR, sRI, sI;
                sR = s * R;
                sRI = s * R * I;
                sI = s * I;
                for(int r=0; r<R; r++){
                    current_rho_immunity_failure[nSR + sR + r] =
                        rho_immunity_failure_values[nSIR + sRI + iR + r];
                }
                current_sigma_immunity_failure[nS + s] =
                    sigma_immunity_failure_values[nSI + sI + i];
            }
            requesting_immunity_update[n] = 0;
        }
    }
    return;
}

void transmission
(
    int S, // number_of_strains
    int L, // number_of_locations
    int N, // number_of_agents
    int id,
    double facemask_transmission_multiplier,
    double current_region_transmission_multiplier,
    const int * current_strain,
    const double * current_disease,
    const int * current_facemask,
    const int * current_location,
    const int * current_region,
    const double * current_infectiousness,
    const double * location_transmission_multiplier,
    const double * mutation_matrix,
    const double * current_sigma_immunity_failure,
    int * infection_event,
    uint64_t * random_state,
    int * random_p
)
{
    double * transmission_force = (double *)malloc(sizeof(double) * L);
    double * sum = (double *)malloc(sizeof(double) * L);
    double * sum_by_strain = (double *)malloc(sizeof(double) * (L * S));
    for(int m=0; m<L; m++){
        transmission_force[m] = 1.0;
        sum[m] = 0.0;
        for(int s=0; s<S; s++){sum_by_strain[(m * S) + s] = 0.0;}
    }

    // Transmission out
    for(int n=0; n<N; n++){
        if(current_region[n] == id && current_strain[n] != -1){
            double f;
            f = current_region_transmission_multiplier *
                location_transmission_multiplier[current_location[n]] *
                (1 + (current_facemask[n] * (facemask_transmission_multiplier - 1))) *
                current_infectiousness[n];
            sum_by_strain[(current_location[n] * S) + current_strain[n]] += f;
            transmission_force[current_location[n]] *= 1 - f;
        }
    }
    for(int m=0; m<L; m++){
        transmission_force[m] = 1 - transmission_force[m];
        for(int s=0; s<S; s++){
            sum[m] += sum_by_strain[(m * S) + s];
        }
    }

    // Transmission in
    double * weights = (double *)malloc(sizeof(double) * S);
    double sum_of_weights;
    int s1, s2;
    for(int n=0; n<N; n++){
        if(current_region[n] == id && current_strain[n] == -1 && current_disease[n] < 1.0){
            double facemask_multiplier;
            facemask_multiplier =
                1 + (current_facemask[n] * (facemask_transmission_multiplier - 1));
            double prob;
            prob = facemask_multiplier * transmission_force[current_location[n]];
            if(bernoulli(random_state, random_p, prob) == 1){
                sum_of_weights = 0;
                for(int s=0; s<S; s++){
                    weights[s] = sum_by_strain[(current_location[n] * S) + s];
                    sum_of_weights += weights[s];
                }
                s1 = random_choice(random_state, random_p, weights, sum_of_weights);
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
        }
    }
    free(transmission_force);
    free(sum);
    free(sum_by_strain);
    free(weights);
    return;
}

void infect
(
    int t,
    int N, // number_of_agents
    int R, // number_of_rho_immunity_outcomes
    int S, // number_of_strains
    int H, // max_preset_length_health
    int G, // max_preset_length_immunity
    int I, // immunity_length
    int immunity_period_ticks,
    const double * current_rho_immunity_failure,
    int * infection_event,
    int * infectiousness_partitions,
    double * infectiousness_values,
    int * infectiousness_lengths,
    int * infectiousness_indexes,
    int * disease_partitions,
    double * disease_values,
    int * disease_lengths,
    int * disease_indexes,
    int * strain_partitions,
    int * strain_values,
    int * strain_lengths,
    int * strain_indexes,
    double * rho_immunity_failure_values,
    double * sigma_immunity_failure_values,
    int * requesting_immunity_update,
    const int * preset_infectiousness_partitions,
    const double * preset_infectiousness_values,
    const int * preset_infectiousness_lengths,
    const int * preset_disease_partitions,
    const double * preset_disease_values,
    const int * preset_disease_lengths,
    const int * preset_strain_partitions,
    const int * preset_strain_values,
    const int * preset_strain_lengths,
    const int * preset_rho_immunity_failure_partitions,
    const double * preset_rho_immunity_failure_values,
    const int * preset_rho_immunity_failure_lengths,
    const int * preset_sigma_immunity_failure_partitions,
    const double * preset_sigma_immunity_failure_values,
    const int * preset_sigma_immunity_failure_lengths,
    const int * presets,
    uint64_t * random_state,
    int * random_p
)
{
    for(int n=0; n<N; n++){
        if(infection_event[n] != -1){

            // The strain with which agent n is to be infected
            int s1;
            s1 = infection_event[n];

            // Reset counters
            infectiousness_indexes[n] = 1;
            disease_indexes[n] = 1;
            strain_indexes[n] = 1;

            // Determine the rho immunity outcome of agent n
            int q, r1, r2;
            q = presets[n];
            r2 = 0;
            while(bernoulli(random_state, random_p,
                    current_rho_immunity_failure[(n * S * R) + (s1 * R) + r2]) == 1){
                r2 += 1;
            }
            r1 = r2;

            // Determine new infectiousness

            infectiousness_lengths[n] =
                preset_infectiousness_lengths[(q * R * S) +
                                                (r1 * R) +
                                                s1];

            for(int i=0; i<infectiousness_lengths[n]; i++){

                infectiousness_values[(n * H) + i] =
                    preset_infectiousness_values[(q * R * S * H) +
                                                    (r1 * S * H) +
                                                    (s1 * H) +
                                                    i];

                infectiousness_partitions[(n * H) + i] =
                    preset_infectiousness_partitions[(q * R * S * H) +
                                                        (r1 * S * H) +
                                                        (s1 * H) +
                                                        i] + t;

            }

            infectiousness_partitions[(n * H) + 0] = -1;

            // Determine new disease

            disease_lengths[n] =
                preset_disease_lengths[(q * R * S) +
                                        (r1 * R) +
                                        s1];

            for(int i=0; i<disease_lengths[n]; i++){

                disease_values[(n * H) + i] =
                    preset_disease_values[(q * R * S * H) +
                                            (r1 * S * H) +
                                            (s1 * H) +
                                            i];

                disease_partitions[(n * H) + i] =
                    preset_disease_partitions[(q * R * S * H) +
                                                (r1 * S * H) +
                                                (s1 * H) +
                                                i] + t;

            }

            disease_partitions[(n * H) + 0] = -1;

            // Determine new strain

            strain_lengths[n] =
                preset_strain_lengths[(q * R * S) +
                                        (r1 * R) +
                                        s1];

            for(int i=0; i<strain_lengths[n]; i++){

                strain_values[(n * H) + i] =
                    preset_strain_values[(q * R * S * H) +
                                            (r1 * S * H) +
                                            (s1 * H) +
                                            i];

                strain_partitions[(n * H) + i] =
                    preset_strain_partitions[(q * R * S * H) +
                                                (r1 * S * H) +
                                                (s1 * H) +
                                                i] + t;
            }

            strain_partitions[(n * H) + 0] = -1;

            // Determine new rho immunity

            int * rho_part = (int *)malloc(sizeof(int) * G);
            double * rho_values = (double *)malloc(sizeof(double) * G * R);
            int rho_length;

            for(int s2=0; s2<S; s2++){

                for(int i=0; i<G; i++){
                    for(int r=0; r<R; r++){

                        rho_values[(i * R) + r] =
                            preset_rho_immunity_failure_values[(q * R * S * S * G * R) +
                                                                (r1 * S * S * G * R) +
                                                                (s1 * S * G * R) +
                                                                (s2 * G * R) +
                                                                (i * R) +
                                                                r];

                    }

                    rho_part[i] =
                        preset_rho_immunity_failure_partitions[(q * R * S * S * G) +
                                                                (r1 * S * S * G) +
                                                                (s1 * S * G) +
                                                                (s2 * G) +
                                                                i] + t;
                    rho_part[0] = -1;

                }

                rho_length =
                    preset_rho_immunity_failure_lengths[(q * R * S * S) +
                                                        (r1 * S * S) +
                                                        (s1 * S) +
                                                        s2];

                for(int i=0; i<I; i++){
                    for(int r=0; r<R; r++){

                        rho_immunity_failure_values[(n * S * I * R) +
                                                    (s2 * I * R) +
                                                    (i * R) +
                                                    r] =

                            rho_immunity_failure_values[(n * S * I * R) +
                                                        (s2 * I * R) +
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

            int * sigma_part = (int *)malloc(sizeof(int) * G);
            double * sigma_values = (double *)malloc(sizeof(double) * G);
            int sigma_length;

            for(int s2=0; s2<S; s2++){

                for(int i=0; i<G; i++){

                    sigma_values[i] =
                        preset_sigma_immunity_failure_values[(q * R * S * S * G) +
                                                                (r1 * S * S * G) +
                                                                (s1 * S * G) +
                                                                (s2 * G) +
                                                                i];

                    sigma_part[i] =
                        preset_sigma_immunity_failure_partitions[(q * R * S * S * G) +
                                                                    (r1 * S * S * G) +
                                                                    (s1 * S * G) +
                                                                    (s2 * G) +
                                                                    i] + t;
                    sigma_part[0] = -1;

                }

                sigma_length =
                    preset_sigma_immunity_failure_lengths[(q * R * S * S) +
                                                            (r1 * S * S) +
                                                            (s1 * S) +
                                                            s2];

                for(int i=0; i<I; i++){

                    sigma_immunity_failure_values[(n * S * I) +
                                                    (s2 * I) +
                                                    i] =

                        sigma_immunity_failure_values[(n * S * I) +
                                                        (s2 * I) +
                                                        i] *

                        sigma_evaluate(&sigma_part[0], &sigma_values[0], sigma_length,
                                        (i + 1) * immunity_period_ticks);

                }
            }
            free(sigma_part);
            free(sigma_values);

            // Complete infection event
            infection_event[n] = -1;
            requesting_immunity_update[n] = 1;
        }
    }
    return;
}
