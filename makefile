CC :=gcc
CFLAGS := -fPIC -shared


all:
	$(CC) $(CFLAGS) -o ./src/pandemia/components/travel_model/default_travel_model_functions  ./src/pandemia/components/travel_model/default_travel_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/components/testing_and_contact_tracing_model/default_testing_and_contact_tracing_model_functions ./src/pandemia/components/testing_and_contact_tracing_model/default_testing_and_contact_tracing_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/components/vaccination_model/default_vaccination_model_functions ./src/pandemia/components/vaccination_model/default_vaccination_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/components/hospitalization_and_death_model/default_hospitalization_and_death_model_functions ./src/pandemia/components/hospitalization_and_death_model/default_hospitalization_and_death_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/components/movement_model/default_movement_model_functions ./src/pandemia/components/movement_model/default_movement_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/components/health_model/default_health_model_functions ./src/pandemia/components/health_model/default_health_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/simulator_functions ./src/pandemia/simulator_functions.c

clean: 
	rm ./src/pandemia/components/travel_model/default_travel_model_functions
	rm ./src/pandemia/components/testing_and_contact_tracing_model/default_testing_and_contact_tracing_model_functions
	rm ./src/pandemia/components/vaccination_model/default_vaccination_model_functions
	rm ./src/pandemia/components/hospitalization_and_death_model/default_hospitalization_and_death_model_functions
	rm ./src/pandemia/components/movement_model/default_movement_model_functions
	rm ./src/pandemia/components/health_model/default_health_model_functions
	rm ./src/pandemia/simulator_functions