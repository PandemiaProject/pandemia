CC :=gcc
CFLAGS := -fPIC -shared

ifeq ($(OS),Windows_NT)
	EXT := .dll
else
	EXT := .so
endif
all:
	$(CC) $(CFLAGS) -o ./src/pandemia/components/travel_model/default_travel_model_functions$(EXT)  ./src/pandemia/components/travel_model/default_travel_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/components/testing_and_contact_tracing_model/default_testing_and_contact_tracing_model_functions$(EXT) ./src/pandemia/components/testing_and_contact_tracing_model/default_testing_and_contact_tracing_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/components/vaccination_model/default_vaccination_model_functions$(EXT) ./src/pandemia/components/vaccination_model/default_vaccination_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/components/hospitalization_and_death_model/default_hospitalization_and_death_model_functions$(EXT) ./src/pandemia/components/hospitalization_and_death_model/default_hospitalization_and_death_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/components/movement_model/default_movement_model_functions$(EXT) ./src/pandemia/components/movement_model/default_movement_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/components/health_model/default_health_model_functions$(EXT) ./src/pandemia/components/health_model/default_health_model_functions.c
	$(CC) $(CFLAGS) -o ./src/pandemia/simulator_functions$(EXT) ./src/pandemia/simulator_functions.c

clean: 
	rm ./src/pandemia/components/travel_model/default_travel_model_functions$(EXT) 
	rm ./src/pandemia/components/testing_and_contact_tracing_model/default_testing_and_contact_tracing_model_functions$(EXT) 
	rm ./src/pandemia/components/vaccination_model/default_vaccination_model_functions$(EXT) 
	rm ./src/pandemia/components/hospitalization_and_death_model/default_hospitalization_and_death_model_functions$(EXT) 
	rm ./src/pandemia/components/movement_model/default_movement_model_functions$(EXT) 
	rm ./src/pandemia/components/health_model/default_health_model_functions$(EXT) 
	rm ./src/pandemia/simulator_functions$(EXT) 