#include <stdlib.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>

#define int32_t_max 2147483647

static inline uint32_t rotl(const uint32_t x, int32_t k) {
	return (x << k) | (x >> (32 - k));
}

// https://prng.di.unimi.it/xoshiro128plusplus.c
uint32_t next(uint32_t * s) {
	const uint32_t result = rotl(s[0] + s[3], 7) + s[0];
	const uint32_t t = s[1] << 9;
	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];
	s[2] ^= t;
	s[3] = rotl(s[3], 11);
	return result;
}

int32_t randrange(uint32_t * s, int32_t num) {
    // Returns an int in the range [0, num)
    return (int32_t) (next(s) % num);
}

int32_t bernoulli(uint32_t * s, double prob_success) {
    // Returns the value 1 if success, 0 otherwise
    int32_t rnd = randrange(s, int32_t_max);
    int32_t result = 0;
    if(rnd < (int32_t) (prob_success * int32_t_max)){
        result = 1;
    }
    return result;
}

int32_t random_choice(uint32_t * s, const double * weights, double sum_of_weights) {
    // Returns an int in the range [0,l) where l is the length of weights
    int32_t rnd = randrange(s, int32_t_max);
    int32_t index = 0;
    while (rnd >= (int32_t) ((weights[index] / sum_of_weights) * int32_t_max)){
        rnd = rnd - (int32_t) (weights[index] / sum_of_weights * int32_t_max);
        ++index;
    }
    return index;
}

void random_shuffle(uint32_t * s, int32_t * array, int32_t n) {
    int32_t i, j, tmp;
    for (i = n - 1; i > 0; i--){
        j = randrange(s, i + 1);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
}

void random_sample(uint32_t * s, int32_t * sample, int32_t n, int32_t * population, int32_t N) {
    int32_t t = 0; // total input records dealt with
    int32_t m = 0; // number of items selected so far
    int32_t u;
    while (m < n)
    {
        u = randrange(s, N - t);
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

double rho_evaluate(const int32_t * partition, const double * values, int32_t length, int32_t r, int32_t R, int32_t t){
    // Evaluates a rho immunity preset on rho immunity outcome r at time t, where R denotes the
    // total number of rho immunity outcomes
    double result;
    if(t < partition[1]){
        result = values[(0 * R) + r];
    } else {
        int32_t j;
        j = length - 1;
        while(t < partition[j]){
            j = j - 1;
        }
        result = values[(j * R) + r];
    }
    return result;
}

double sigma_evaluate(const int32_t * partition, const double * values, int32_t length, int32_t t){
    // Evaluates a sigma immunity preset at time t
    double result;
    if(t < partition[1]){
        result = values[0];
    } else {
        int32_t j;
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
                int32_t t,
                int32_t N, // number_of_agents
                int32_t S, // number_of_strains
                int32_t I, // immunity_length
                int32_t R, // number_of_rho_immunity_outcomes
                int32_t H, // max_preset_length_health
                int32_t immunity_period_ticks,
                const int32_t * i_partitions,
                int32_t * i_indexes,
                const double * i_values,
                const int32_t * i_lengths,
                double * i_current,
                const int32_t * d_partitions,
                int32_t * d_indexes,
                const double * d_values,
                const int32_t * d_lengths,
                double * d_current,
                const int32_t * s_partitions,
                int32_t * s_indexes,
                const int32_t * s_values,
                const int32_t * s_lengths,
                int32_t * s_current,
                const double * rho_immunity_failure_values,
                double * current_rho_immunity_failure,
                const double * sigma_immunity_failure_values,
                double * current_sigma_immunity_failure,
                int32_t * requesting_immunity_update
        )
{
    if(t % immunity_period_ticks == 0){for(int32_t n=0; n<N; n++){requesting_immunity_update[n] = 1;}}
    int32_t i = t / immunity_period_ticks;
    for(int32_t n=0; n<N; n++){
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
            int32_t nSR, nSIR, nS, nSI, iR;
            iR = i * R;
            nSR = n * S * R;
            nSIR = n * S * I * R;
            nS = n * S;
            nSI = n * S * I;
            for(int32_t s=0; s<S; s++){
                int32_t sR, sRI, sI;
                sR = s * R;
                sRI = s * R * I;
                sI = s * I;
                for(int32_t r=0; r<R; r++){
                    current_rho_immunity_failure[nSR + sR + r] =
                            rho_immunity_failure_values[nSIR + sRI + iR + r];
                }
                current_sigma_immunity_failure[nS + s] =
                        sigma_immunity_failure_values[nSI + sI + i];
            }
            requesting_immunity_update[n] = 0;
        }
    }
}

void transmission
        (
                int32_t S, // number_of_strains
                int32_t L, // number_of_locations
                int32_t N, // number_of_agents
                int32_t id,
                int32_t sir_rescaling_int,
                int32_t ticks_in_day,
                int32_t A, // number of age groups
                const double * beta,
                const int32_t * subpopulation_index,
                const double * subpopulation_mixing_matrix,
                double facemask_transmission_multiplier,
                double current_region_transmission_multiplier,
                const int32_t * current_strain,
                const double * current_disease,
                const int32_t * current_facemask,
                const int32_t * current_location,
                const int32_t * current_region,
                const double * current_infectiousness,
                const double * location_transmission_multiplier,
                const double * mutation_matrix,
                const double * current_sigma_immunity_failure,
                int32_t * infection_event,
                uint32_t * random_state
        )
{
    double * transmission_force_by_age_group = (double *)malloc(sizeof(double) * L * A);
    double * sum_by_strain_by_age_group = (double *)malloc(sizeof(double) * L * S * A);
    int32_t * num_agents_by_location_by_age_group = (int32_t *)malloc(sizeof(int32_t) * L * A);
    for(int32_t m=0; m<L; m++){
        for(int32_t a=0; a<A; a++){
            transmission_force_by_age_group[(m * A) + a] = 1.0;
            num_agents_by_location_by_age_group[(m * A) + a] = 0;
            for(int32_t s=0; s<S; s++){sum_by_strain_by_age_group[(m * S * A) + (s * A) + a] = 0.0;}
        }
    }

    // Transmission out
    if(sir_rescaling_int == 1){
        for(int32_t n=0; n<N; n++){
            if(current_region[n] == id){
                num_agents_by_location_by_age_group[(current_location[n] * A) +
                                                    subpopulation_index[n]] += 1;
            }
        }
    }
    for(int32_t a=0; a<A; a++){
        for(int32_t n=0; n<N; n++){
            if(current_region[n] == id && current_strain[n] != -1){
                double f;
                f = current_region_transmission_multiplier *
                    location_transmission_multiplier[current_location[n]] *
                    (1 + (current_facemask[n] * (facemask_transmission_multiplier - 1))) *
                    current_infectiousness[n] * beta[current_strain[n]];
                if(sir_rescaling_int == 1){
                    f *= subpopulation_mixing_matrix[(a * A) + subpopulation_index[n]] /
                         (num_agents_by_location_by_age_group[(current_location[n] * A) +
                                                              subpopulation_index[n]] * ticks_in_day);
                }
                sum_by_strain_by_age_group[(current_location[n] * S * A) +
                                           (current_strain[n] * A) + a] += f;
                transmission_force_by_age_group[(current_location[n] * A) + a] *= 1 - f;
            }
        }
    }
    for(int32_t a=0; a<A; a++){
        for(int32_t m=0; m<L; m++){
            transmission_force_by_age_group[(m * A) + a] =
                    1 - transmission_force_by_age_group[(m * A) + a];
        }
    }

    // Transmission in
    double * weights = (double *)malloc(sizeof(double) * S);
    double sum_of_weights;
    int32_t s1, s2;
    for(int32_t n=0; n<N; n++){
        if(current_region[n] == id && current_strain[n] == -1 && current_disease[n] < 1.0){
            double facemask_multiplier;
            facemask_multiplier =
                    1 + (current_facemask[n] * (facemask_transmission_multiplier - 1));
            double prob;
            if(S > 1){
                prob = facemask_multiplier *
                       transmission_force_by_age_group[(current_location[n] * A) +
                                                       subpopulation_index[n]];
                if(bernoulli(random_state, prob) == 1){
                    sum_of_weights = 0;
                    for(int32_t s=0; s<S; s++){
                        weights[s] = sum_by_strain_by_age_group[(current_location[n] * S * A) +
                                                                (s * A) + subpopulation_index[n]];
                        sum_of_weights += weights[s];
                    }
                    s1 = random_choice(random_state, weights, sum_of_weights);
                    prob = current_sigma_immunity_failure[(n * S) + s1];
                    if(bernoulli(random_state, prob) == 1){
                        sum_of_weights = 0;
                        for(int32_t s=0; s<S; s++){
                            weights[s] = mutation_matrix[(s1 * S) + s];
                            sum_of_weights += weights[s];
                        }
                        s2 = random_choice(random_state, weights, sum_of_weights);
                        infection_event[n] = s2;
                    }
                }
            } else {
                s1 = 0;
                prob = facemask_multiplier *
                       transmission_force_by_age_group[(current_location[n] * A) +
                                                       subpopulation_index[n]] *
                       current_sigma_immunity_failure[(n * S) + s1];
                if(bernoulli(random_state, prob) == 1){
                    infection_event[n] = s1;
                }
            }
        }
    }
    free(transmission_force_by_age_group);
    free(sum_by_strain_by_age_group);
    free(num_agents_by_location_by_age_group);
    free(weights);
}

void infect
        (
                int32_t t,
                int32_t N, // number_of_agents
                int32_t R, // number_of_rho_immunity_outcomes
                int32_t S, // number_of_strains
                int32_t H, // max_preset_length_health
                int32_t G, // max_preset_length_immunity
                int32_t I, // immunity_length
                int32_t immunity_period_ticks,
                const double * current_rho_immunity_failure,
                int32_t * infection_event,
                int32_t * infectiousness_partitions,
                double * infectiousness_values,
                int32_t * infectiousness_lengths,
                int32_t * infectiousness_indexes,
                int32_t * disease_partitions,
                double * disease_values,
                int32_t * disease_lengths,
                int32_t * disease_indexes,
                int32_t * strain_partitions,
                int32_t * strain_values,
                int32_t * strain_lengths,
                int32_t * strain_indexes,
                double * rho_immunity_failure_values,
                double * sigma_immunity_failure_values,
                int32_t * requesting_immunity_update,
                const int32_t * preset_infectiousness_partitions,
                const double * preset_infectiousness_values,
                const int32_t * preset_infectiousness_lengths,
                const int32_t * preset_disease_partitions,
                const double * preset_disease_values,
                const int32_t * preset_disease_lengths,
                const int32_t * preset_strain_partitions,
                const int32_t * preset_strain_values,
                const int32_t * preset_strain_lengths,
                const int32_t * preset_rho_immunity_failure_partitions,
                const double * preset_rho_immunity_failure_values,
                const int32_t * preset_rho_immunity_failure_lengths,
                const int32_t * preset_sigma_immunity_failure_partitions,
                const double * preset_sigma_immunity_failure_values,
                const int32_t * preset_sigma_immunity_failure_lengths,
                const int32_t * presets,
                uint32_t * random_state
        )
{
    for(int32_t n=0; n<N; n++){
        if(infection_event[n] != -1){

            // The strain with which agent n is to be infected
            int32_t s1;
            s1 = infection_event[n];

            // Reset counters
            infectiousness_indexes[n] = 1;
            disease_indexes[n] = 1;
            strain_indexes[n] = 1;

            // Determine the rho immunity outcome of agent n
            int32_t q, r1, r2;
            q = presets[n];
            r2 = 0;
            while(bernoulli(random_state,
                            current_rho_immunity_failure[(n * S * R) + (s1 * R) + r2]) == 1){
                r2 += 1;
            }
            r1 = r2;
            // Determine new infectiousness

            infectiousness_lengths[n] =
                    preset_infectiousness_lengths[(q * R * S) +
                                                  (r1 * S) +
                                                  s1];

            for(int32_t i=0; i<infectiousness_lengths[n]; i++){

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
                                           (r1 * S) +
                                           s1];

            for(int32_t i=0; i<disease_lengths[n]; i++){

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
                                          (r1 * S) +
                                          s1];

            for(int32_t i=0; i<strain_lengths[n]; i++){

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

            int32_t * rho_part = (int32_t *)malloc(sizeof(int32_t) * G);
            double * rho_values = (double *)malloc(sizeof(double) * G * R);
            int32_t rho_length;

            for(int32_t s2=0; s2<S; s2++){

                for(int32_t i=0; i<G; i++){
                    for(int32_t r=0; r<R; r++){

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

                for(int32_t i=0; i<I; i++){
                    for(int32_t r=0; r<R; r++){

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

            int32_t * sigma_part = (int32_t *)malloc(sizeof(int32_t) * G);
            double * sigma_values = (double *)malloc(sizeof(double) * G);
            int32_t sigma_length;

            for(int32_t s2=0; s2<S; s2++){

                for(int32_t i=0; i<G; i++){

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

                for(int32_t i=0; i<I; i++){

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
}

void dynamics_hospitalization_and_death
        (
                int32_t N, // number_of_agents
                int32_t number_of_hospitals,
                int32_t number_of_cemeteries,
                int32_t id,
                double hospital_threshold,
                const int32_t * current_location,
                const int32_t * current_region,
                const double * current_disease,
                const int32_t * requesting_location_update,
                int32_t * requested_location_update,
                int32_t * requesting_facemask_update,
                int32_t * requested_facemask_update,
                int32_t * in_hospital,
                const int32_t * hospitals,
                int32_t * in_cemetery,
                const int32_t * cemeteries,
                uint32_t * random_state
        )
{
    // Hospitals
    if(number_of_hospitals > 0){
        for(int32_t n=0; n<N; n++){
            if(current_region[n] == id && requesting_location_update[n] == 1){
                if(current_disease[n] > hospital_threshold){
                    if(in_hospital[n] == 0){
                        int32_t index = 0;
                        index = randrange(random_state, number_of_hospitals);
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
        for(int32_t n=0; n<N; n++){
            if(current_region[n] == id && requesting_location_update[n] == 1){
                if(current_disease[n] == 1.0){
                    if(in_cemetery[n] == 0){
                        int32_t index = 0;
                        index = randrange(random_state, number_of_cemeteries);
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


void update_movement
        (
                int32_t N, // number_of_agents
                int32_t * requesting_location_update,
                const int32_t * requested_location_update,
                int32_t * current_location,
                int32_t * requesting_facemask_update,
                const int32_t * requested_facemask_update,
                int32_t * current_facemask
        )
{
    for(int32_t n=0; n<N; n++){
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
                int32_t N, // number_of_agents
                int32_t A, // number_of_activities
                int32_t lockdown_intervention,
                int32_t facemask_intervention,
                int32_t id,
                int32_t use_weights,
                int32_t t,
                int32_t offset,
                int32_t ticks_in_week,
                int32_t max_num_activity_locations,
                const int32_t * current_region,
                const uint8_t * weekly_routines,
                const int32_t * current_facemask,
                const int32_t * wears_facemask,
                int32_t * requested_facemask_update,
                int32_t * requesting_facemask_update,
                const int32_t * activity_locations,
                const double * activity_location_weights,
                const int32_t * num_activity_locations,
                const int32_t * location_closure,
                int32_t * requested_location_update,
                const int32_t * home_location,
                const int32_t * current_quarantine,
                int32_t * requesting_location_update,
                uint32_t * random_state
        )
{
    int32_t t_now = (t + offset) % ticks_in_week;
    int32_t t_next = (t + 1 + offset) % ticks_in_week;

    // Determine new locations, for those agents starting a new activity
    for(int32_t n=0; n<N; n++){
        if(current_region[n] == id){
            uint8_t new_activity = weekly_routines[(n * ticks_in_week) + t_next];
            if(new_activity != weekly_routines[(n * ticks_in_week) + t_now]){
                // If the agent is not currently wearing a facemask but should be, put it on
                if(current_facemask[n] == 0 && facemask_intervention == 1 &&
                   wears_facemask[(n * A) + new_activity] == 1){
                    requested_facemask_update[n] = 1;
                    requesting_facemask_update[n] = 1;
                }
                // If the agent is currently wearing a facemask but shouldn't be, take it off
                if(current_facemask[n] == 1 && (facemask_intervention == 0 ||
                                                wears_facemask[(n * A) + new_activity] == 0)){
                    requested_facemask_update[n] = 0;
                    requesting_facemask_update[n] = 1;
                }
                // Randomly select an index to determine the new location
                int32_t index;
                int32_t num = num_activity_locations[(n * A) + new_activity];
                if(use_weights == 1){
                    double * weights = (double *)malloc(sizeof(double) * num);
                    double sum_of_weights = 0;
                    int32_t index_0 = (n * A * max_num_activity_locations) +
                                  (new_activity * max_num_activity_locations);
                    for(int32_t j=0; j<num; j++){
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

int32_t test
        (
                uint32_t * s,
                double infectiousness,
                double test_threshold,
                double test_false_negative
        )
{
    // Tests agent based on infectiousness
    int32_t test_result = 0;
    if(infectiousness >= test_threshold){
        test_result = bernoulli(s, 1 - test_false_negative);
    }
    return test_result;
}

void default_testing_and_contact_tracing_dynamics
        (
                int32_t N, // number_of_agents
                int32_t day,
                int32_t quarantine_period_days,
                int32_t num_to_test_random,
                int32_t num_to_test_symptomatic,
                int32_t num_to_test_contact,
                int32_t max_regular_contacts_to_test,
                int32_t id,
                double symptomatic_disease_threshold,
                double test_threshold,
                double test_false_negative,
                double prob_quarantine_with_symptoms_without_test,
                double prob_quarantine_with_contact_without_test,
                const int32_t * num_regular_contacts_to_test,
                const int32_t * regular_contacts_to_test,
                const int32_t * current_region,
                const double * current_infectiousness,
                const double * current_disease,
                double * yesterdays_disease,
                int32_t * end_of_quarantine_days,
                int32_t * current_quarantine,
                uint32_t * random_state
        )
{
    // End quarantine if necessary
    for(int32_t n=0; n<N; n++){
        if(end_of_quarantine_days[n] == day){
            current_quarantine[n] = 0;
        }
    }

    // Population eligible to be tested today
    int32_t * eligible = (int32_t *)malloc(sizeof(int32_t) * N);
    for(int32_t n=0; n<N; n++){
        if((current_region[n] == id) && (current_disease[n] < 1.0)){
            eligible[n] = 1;
        } else {
            eligible[n] = 0;
        }
    }

    // A record of agents newly testing positive today
    int32_t * newly_testing_positive = (int32_t *)malloc(sizeof(int32_t) * N);
    for(int32_t n=0; n<N; n++){newly_testing_positive[n] = 0;}

    // Random testing (test eligible agents at random)
    if(num_to_test_random > 0){
        int32_t num_eligible_rand = 0;
        for(int32_t n=0; n<N; n++){if(eligible[n] == 1){num_eligible_rand += 1;}}
        int32_t * eligible_agents_rand = (int32_t *)malloc(sizeof(int32_t) * num_eligible_rand);
        int32_t j = 0;
        for(int32_t n=0; n<N; n++){if(eligible[n] == 1){eligible_agents_rand[j] = n; j += 1;}}
        int32_t num_agents_to_test_random;
        num_agents_to_test_random = fmin(num_to_test_random, num_eligible_rand);
        int32_t * agents_to_test_random = (int32_t *)malloc(sizeof(int32_t) * num_agents_to_test_random);
        random_sample(random_state, agents_to_test_random, num_agents_to_test_random,
                      eligible_agents_rand, num_eligible_rand);
        for(int32_t j=0; j<num_agents_to_test_random; j++){
            int32_t n, test_result;
            n = agents_to_test_random[j];
            eligible[n] = 0;
            test_result = test(random_state, current_infectiousness[n], test_threshold,
                               test_false_negative);
            if(test_result == 1){
                newly_testing_positive[n] = 1;
                current_quarantine[n] = 1;
                end_of_quarantine_days[n] = day + quarantine_period_days;
            }
        }
        free(agents_to_test_random);
        free(eligible_agents_rand);
    }

    // Symptomatic testing (test eligible agents who have just become symptomatic)
    int32_t num_eligible_symp = 0;
    for(int32_t n=0; n<N; n++){if(eligible[n] == 1){num_eligible_symp += 1;}}
    int32_t * eligible_agents_symp = (int32_t *)malloc(sizeof(int32_t) * num_eligible_symp);
    int32_t j = 0;
    for(int32_t n=0; n<N; n++){if(eligible[n] == 1){eligible_agents_symp[j] = n; j += 1;}}
    random_shuffle(random_state, eligible_agents_symp, num_eligible_symp);
    for(int32_t j=0; j<num_eligible_symp; j++){
        int32_t n, test_result;
        n = eligible_agents_symp[j];
        if(current_disease[n] >= symptomatic_disease_threshold){
            if(yesterdays_disease[n] < symptomatic_disease_threshold){
                if(num_to_test_symptomatic > 0){
                    test_result = test(random_state, current_infectiousness[n], test_threshold,
                                       test_false_negative);
                    if(test_result == 1){
                        newly_testing_positive[n] = 1;
                        current_quarantine[n] = 1;
                        end_of_quarantine_days[n] = day + quarantine_period_days;
                    }
                    num_to_test_symptomatic -= 1;
                } else {
                    if(bernoulli(random_state, prob_quarantine_with_symptoms_without_test)){
                        current_quarantine[n] = 1;
                        end_of_quarantine_days[n] = day + quarantine_period_days;
                    }
                }
            }
        }
    }
    free(eligible_agents_symp);

    // Contact tracing (contact trace and test eligible agents)
    int32_t num_newly_testing_positive = 0;
    for(int32_t n=0; n<N; n++){if(newly_testing_positive[n] == 1){num_newly_testing_positive += 1;}}
    if((num_newly_testing_positive > 0) && (num_to_test_contact > 0)){
        int32_t * at_risk = (int32_t *)malloc(sizeof(int32_t) * N);
        for(int32_t n=0; n<N; n++){at_risk[n] = 0;}
        for(int32_t n1=0; n1<N; n1++){
            if(newly_testing_positive[n1] == 1){
                for(int32_t j=0; j<num_regular_contacts_to_test[n1]; j++){
                    int32_t n2;
                    n2 = regular_contacts_to_test[(n1 * max_regular_contacts_to_test) + j];
                    if((current_region[n2] == id) && (current_disease[n2] < 1.0)){
                        at_risk[n2] = 1;
                    }
                }
            }
        }
        for(int32_t n=0; n<N; n++){if(newly_testing_positive[n] == 1){at_risk[n] = 0;}}
        int32_t num_agents_at_risk = 0;
        for(int32_t n=0; n<N; n++){if(at_risk[n] == 1){num_agents_at_risk += 1;}}
        int32_t * agents_at_risk = (int32_t *)malloc(sizeof(int32_t) * num_agents_at_risk);
        int32_t j = 0;
        for(int32_t n=0; n<N; n++){if(at_risk[n] == 1){agents_at_risk[j] = n; j += 1;}}
        random_shuffle(random_state, agents_at_risk, num_agents_at_risk);
        for(int32_t j=0; j<fmin(num_to_test_contact, num_agents_at_risk); j++){
            int32_t n, test_result;
            n = agents_at_risk[j];
            test_result = test(random_state, current_infectiousness[n], test_threshold,
                               test_false_negative);
            if(test_result == 1){
                current_quarantine[n] = 1;
                end_of_quarantine_days[n] = day + quarantine_period_days;
            }
        }
        for(int32_t j=fmin(num_to_test_contact, num_agents_at_risk); j<num_agents_at_risk; j++){
            int32_t n;
            n = agents_at_risk[j];
            if(bernoulli(random_state, prob_quarantine_with_contact_without_test)){
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
    for(int32_t n=0; n<N; n++){
        yesterdays_disease[n] = current_disease[n];
    }
    return;
}

void close_borders
        (
                int32_t R, // number_of_regions
                int32_t id,
                double scale_factor,
                double current_border_closure_multiplier,
                int32_t * agents_travelling_matrix,
                const double * baseline_agents_travelling_matrix
        )
{
    int32_t r1 = id;
    double rescaled;
    for(int32_t r2=0; r2<R; r2++){
        rescaled = baseline_agents_travelling_matrix[(r1 * R) + r2] * scale_factor *
                   current_border_closure_multiplier;
        agents_travelling_matrix[(r1 * R) + r2] = (int32_t) rescaled;
    }
    return;
}

void determine_travellers
        (
                int32_t N, // number_of_agents
                int32_t id,
                int32_t R, // number_of_regions
                const double * current_disease,
                const int32_t * current_strain,
                int32_t * current_region,
                const int32_t * agents_travelling_matrix,
                uint32_t * random_state
        )
{
    int32_t r1 = id;
    int32_t total_num_to_travel = 0;

    // Determine how many agents travel from this region today
    for(int32_t r2=0; r2<R; r2++){total_num_to_travel += agents_travelling_matrix[(r1 * R) + r2];}

    if(total_num_to_travel > 0){

        // Determine who is eligible to travel from this region today and count how many
        int32_t * agents_eligible_to_travel = (int32_t *)malloc(sizeof(int32_t) * N);
        int32_t num_eligible_to_travel = 0;
        for(int32_t n=0; n<N; n++){
            if(current_strain[n] == -1 && current_disease[n] < 1.0){
                agents_eligible_to_travel[num_eligible_to_travel] = n;
                num_eligible_to_travel += 1;
            }
        }

        // Among these eligible agents determine who actually travels today
        int32_t num_agents_to_travel;
        num_agents_to_travel = (int32_t) fmin((double)total_num_to_travel, (double)num_eligible_to_travel);
        int32_t * agents_to_travel = (int32_t *)malloc(sizeof(int32_t) * num_agents_to_travel);
        random_sample(random_state, agents_to_travel, num_agents_to_travel,
                      agents_eligible_to_travel, num_eligible_to_travel);

        // Among those who actually travel determine to which region they travel
        int32_t num_agents_to_travel_r2;
        int32_t j_min, j_max;
        j_min = 0;
        j_max = 0;
        for(int32_t r2=0; r2<R; r2++){
            num_agents_to_travel_r2 = agents_travelling_matrix[(r1 * R) + r2];
            j_max += num_agents_to_travel_r2;
            for(int32_t j=(int32_t)fmin((double)j_min, (double)num_agents_to_travel); j<(int32_t)fmin((double)j_max, (double)num_agents_to_travel); j++){
                current_region[agents_to_travel[j]] = r2;
            }
            j_min += num_agents_to_travel_r2;
        }

        free(agents_to_travel);
        free(agents_eligible_to_travel);
    }
    return;
}

void transmission_out
        (
                int32_t N, // number_of_agents
                int32_t S, // number_of_strains
                int32_t id,
                const double * beta,
                double facemask_transmission_multiplier,
                double travel_multiplier,
                double current_region_transmission_multiplier,
                const int32_t * current_region,
                const double * current_infectiousness,
                const int32_t * current_strain,
                const int32_t * current_facemask,
                double * sum_f_by_strain,
                double * transmission_force
        )
{
    double facemask_multiplier, f;
    for(int32_t n=0; n<N; n++){
        if(current_region[n] == id && current_strain[n] != -1){
            facemask_multiplier = 1 + current_facemask[n] * (facemask_transmission_multiplier - 1);
            f = current_region_transmission_multiplier *
                facemask_multiplier * current_infectiousness[n] * beta[current_strain[n]];
            f = (f * travel_multiplier) / N;
            f = fmin(f, 1.0);
            sum_f_by_strain[(id * S) + current_strain[n]] += f;
            transmission_force[id] *= 1 - f;
        }
    }
    transmission_force[id] = 1 - transmission_force[id];
    return;
}

void transmission_in
        (
                int32_t S, // number_of_strains
                int32_t N, // number_of_agents
                int32_t r1,
                int32_t * current_facemask,
                int32_t * current_region,
                double facemask_transmission_multiplier,
                double * sum_f_by_strain,
                double * current_sigma_immunity_failure,
                int32_t * infection_event,
                double * transmission_force,
                double * mutation_matrix,
                uint32_t * random_state
        )
{
    double * weights = (double *)malloc(sizeof(double) * S);
    double sum_of_weights;
    int32_t s1, s2;
    for(int32_t n=0; n<N; n++){
        if(current_region[n] != r1){
            int32_t r2 = current_region[n];
            double facemask_multiplier;
            facemask_multiplier =
                    1 + current_facemask[n] * (facemask_transmission_multiplier - 1);
            double prob;
            if(S > 1){
                prob = facemask_multiplier * transmission_force[r2];
                if(bernoulli(random_state, prob) == 1){
                    sum_of_weights = 0;
                    for(int32_t s=0; s<S; s++){
                        weights[s] = sum_f_by_strain[(r2 * S) + s];
                        sum_of_weights += weights[s];
                    }
                    s1 = random_choice(random_state, weights, sum_of_weights);
                    double prob;
                    prob = current_sigma_immunity_failure[(n * S) + s1];
                    if(bernoulli(random_state, prob) == 1){
                        sum_of_weights = 0;
                        for(int32_t s=0; s<S; s++){
                            weights[s] = mutation_matrix[(s1 * S) + s];
                            sum_of_weights += weights[s];
                        }
                        s2 = random_choice(random_state, weights, sum_of_weights);
                        infection_event[n] = s2;
                    }
                }
            } else {
                s1 = 0;
                prob = facemask_multiplier *
                       transmission_force[r2] *
                       current_sigma_immunity_failure[(n * S) + s1];
                if(bernoulli(random_state, prob) == 1){
                    infection_event[n] = s1;
                }
            }
        }

    }
    free(weights);
    return;
}

void dynamics_vaccination
        (
                int32_t day,
                int32_t ticks_in_day,
                int32_t N, // number_of_agents
                int32_t S, // number_of_strains
                int32_t I, // immunity_length
                int32_t R, // number_of_rho_immunity_outcomes
                int32_t A, // number_of_age_groups
                int32_t V, // number_of_vaccines
                int32_t W, // max_vaccine_length_immunity
                int32_t booster_waiting_time_days,
                int32_t immunity_period_ticks,
                int32_t id,
                const double * vaccine_rho_immunity_failure_values,
                const int32_t * vaccine_rho_immunity_failure_partitions,
                const int32_t * vaccine_rho_immunity_failure_lengths,
                const double * vaccine_sigma_immunity_failure_values,
                const int32_t * vaccine_sigma_immunity_failure_partitions,
                const int32_t * vaccine_sigma_immunity_failure_lengths,
                const int32_t * num_to_vaccinate,
                double * rho_immunity_failure_values,
                double * sigma_immunity_failure_values,
                const int32_t * current_region,
                const double * current_disease,
                const int32_t * current_strain,
                int32_t * requesting_immunity_update,
                int32_t * most_recent_first_dose,
                int32_t * vaccine_hesitant,
                uint32_t * random_state
        )
{

    for(int32_t age_group_index=0; age_group_index<A; age_group_index++){

        // Determine who is eligible to be vaccinated today for this age group
        int32_t num_eligible = 0;
        int32_t * eligible = (int32_t *)malloc(sizeof(int32_t) * N);
        for(int32_t n=0; n<N; n++){
            if((current_region[n] == id) &&
               (current_disease[n] < 1.0) &&
               (day - most_recent_first_dose[n] >= booster_waiting_time_days) &&
               (current_strain[n] == -1) &&
               (vaccine_hesitant[n] == 0)){
                eligible[n] = 1;
                num_eligible += 1;
            } else {
                eligible[n] = 0;
            }
        }
        int32_t * eligible_agents = (int32_t *)malloc(sizeof(int32_t) * num_eligible);
        int32_t j = 0;
        for(int32_t n=0; n<N; n++){if(eligible[n] == 1){eligible_agents[j] = n; j += 1;}}

        // Vaccinate a sample of the eligible population and update their immunity functions
        int32_t t = day * ticks_in_day;
        int32_t num_agents_to_vaccinate = 0;
        for(int32_t v=0; v<V; v++){
            num_agents_to_vaccinate += num_to_vaccinate[(age_group_index * V) + v];
        }
        num_agents_to_vaccinate = fmin(num_eligible, num_agents_to_vaccinate);
        int32_t * agents_to_vaccinate = (int32_t *)malloc(sizeof(int32_t) * num_agents_to_vaccinate);
        random_sample(random_state, agents_to_vaccinate, num_agents_to_vaccinate, eligible_agents,
                      num_eligible);
        int32_t v = 0;
        int32_t num_vaccinated = 0;
        for(int32_t j=0; j<num_agents_to_vaccinate; j++){
            int32_t n;
            n = agents_to_vaccinate[j];
            if(v < V){

                // Determine new rho immunity

                int32_t * rho_part = (int32_t *)malloc(sizeof(int32_t) * W);
                double * rho_values = (double *)malloc(sizeof(double) * W * R);
                int32_t rho_length;

                for(int32_t s=0; s<S; s++){

                    for(int32_t i=0; i<W; i++){

                        for(int32_t r=0; r<R; r++){

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

                    for(int32_t i=0; i<I; i++){
                        for(int32_t r=0; r<R; r++){

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

                int32_t * sigma_part = (int32_t *)malloc(sizeof(int32_t) * W);
                double * sigma_values = (double *)malloc(sizeof(double) * W);
                int32_t sigma_length;

                for(int32_t s=0; s<S; s++){

                    for(int32_t i=0; i<W; i++){

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
                            vaccine_sigma_immunity_failure_lengths[(v * S) + s];

                    for(int32_t i=0; i<I; i++){

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
            if(num_vaccinated >= num_to_vaccinate[(age_group_index * V) + v]){
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

void safir3_transmission_model
    (
        int32_t L, // number_of_locations
        int32_t N, // number_of_agents
        int32_t id,
        int32_t sir_rescaling_int,
        int32_t ticks_in_day,
        int32_t A, // number of age groups
        const double * beta,
        const int32_t * subpopulation_index,
        const double * subpopulation_mixing_matrix,
        double facemask_transmission_multiplier,
        double current_region_transmission_multiplier,
        const int32_t * current_strain,
        const double * current_disease,
        const int32_t * current_facemask,
        const int32_t * current_location,
        const int32_t * current_region,
        const double * current_infectiousness,
        const double * location_transmission_multiplier,
        const double * ef_transmission,
        const double * ef_infection,
        int32_t * infection_event,
        uint32_t * random_state
    )
{
    double * transmission_force_by_age_group = (double *)malloc(sizeof(double) * L * A);
    int32_t * num_agents_by_location_by_age_group = (int32_t *)malloc(sizeof(int32_t) * L * A);
    for(int32_t m=0; m<L; m++){
        for(int32_t a=0; a<A; a++){
            transmission_force_by_age_group[(m * A) + a] = 1.0;
            num_agents_by_location_by_age_group[(m * A) + a] = 0;
        }
    }

    // Transmission out
    if(sir_rescaling_int == 1){
        for(int32_t n=0; n<N; n++){
            if(current_region[n] == id){
                num_agents_by_location_by_age_group[(current_location[n] * A) +
                                                    subpopulation_index[n]] += 1;
            }
        }
    }
    for(int32_t a=0; a<A; a++){
        for(int32_t n=0; n<N; n++){
            if(current_region[n] == id && current_strain[n] != -1){
                double f;
                f = current_region_transmission_multiplier *
                    location_transmission_multiplier[current_location[n]] *
                    (1 + (current_facemask[n] * (facemask_transmission_multiplier - 1))) *
                    current_infectiousness[n] * beta[0] * ef_transmission[n];
                if(sir_rescaling_int == 1){
                    f *= subpopulation_mixing_matrix[(a * A) + subpopulation_index[n]] /
                         (num_agents_by_location_by_age_group[(current_location[n] * A) +
                                                            subpopulation_index[n]] * ticks_in_day);
                }
                transmission_force_by_age_group[(current_location[n] * A) + a] *= 1 - f;
            }
        }
    }
    for(int32_t a=0; a<A; a++){
        for(int32_t m=0; m<L; m++){
            transmission_force_by_age_group[(m * A) + a] =
                1 - transmission_force_by_age_group[(m * A) + a];
        }
    }

    // Transmission in
    for(int32_t n=0; n<N; n++){
        if(current_region[n] == id && current_strain[n] == -1 && current_disease[n] < 1.0){
            double facemask_multiplier;
            facemask_multiplier =
                1 + (current_facemask[n] * (facemask_transmission_multiplier - 1));
            double prob;
            prob = facemask_multiplier * ef_infection[n] *
                   transmission_force_by_age_group[(current_location[n] * A) +
                                                    subpopulation_index[n]];
            if(bernoulli(random_state, prob) == 1){
                infection_event[n] = 0;
            }
        }
    }
    free(transmission_force_by_age_group);
    free(num_agents_by_location_by_age_group);
    return;
}


void safir_3_transmission_out
    (
        int32_t N, // number_of_agents
        int32_t id,
        const double * beta,
        double facemask_transmission_multiplier,
        double travel_multiplier,
        double current_region_transmission_multiplier,
        const int32_t * current_region,
        const double * current_infectiousness,
        const double * ef_transmission,
        const int32_t * current_strain,
        const int32_t * current_facemask,
        double * transmission_force
    )
{
    double facemask_multiplier, f;
    for(int32_t n=0; n<N; n++){
        if(current_region[n] == id && current_strain[n] != -1){
            facemask_multiplier = 1 + current_facemask[n] * (facemask_transmission_multiplier - 1);
            f = current_region_transmission_multiplier * ef_transmission[n] *
                facemask_multiplier * current_infectiousness[n] * beta[current_strain[n]];
            f = (f * travel_multiplier) / N;
            f = fmin(f, 1.0);
            transmission_force[id] *= 1 - f;
        }
    }
    transmission_force[id] = 1 - transmission_force[id];
    return;
}

void safir3_transmission_in
    (
        int32_t N, // number_of_agents
        int32_t r1,
        int32_t * current_facemask,
        int32_t * current_region,
        double facemask_transmission_multiplier,
        double * ef_infection,
        int32_t * infection_event,
        double * transmission_force,
        uint32_t * random_state
    )
{
    for(int32_t n=0; n<N; n++){
        if(current_region[n] != r1){
            double facemask_multiplier;
            facemask_multiplier = 1 + current_facemask[n] * (facemask_transmission_multiplier - 1);
            double prob;
            prob = facemask_multiplier * transmission_force[current_region[n]] * ef_infection[n];
            if(bernoulli(random_state, prob) == 1){
                infection_event[n] = 0;
            }
        }

    }
    return;
}
