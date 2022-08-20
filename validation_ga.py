
import pygad
import logging
import logging.config
from tqdm import tqdm
import copy
import csv
from multiprocessing import Pool
import matplotlib.pyplot as plt
import numpy as np

from pandemia.utils import instantiate_class
from pandemia.messagebus import MessageBus
from pandemia.sim_factory import SimulationFactory

# Global module log
log = logging.getLogger()

def build_components(sim_factory):
    """The model also features a number of components, representing agent mobility, health and
    number of other features. There is also a scale factor, which can be used to rescale the numbers
    of tests and vaccine doses given as input, the numbers of agents travelling between the regions,
    and the final output."""

    config = sim_factory.config

    scale_factor = config['scale_factor']

    enable_ctypes = config['enable_ctypes']

    # Create seasonal effects model
    seasonal_effects_class = config['seasonal_effects_model.__type__']
    seasonal_effects_config = config.subconfig('seasonal_effects_model')
    seasonal_effects = instantiate_class("pandemia.components.seasonal_effects_model",
                                         seasonal_effects_class, seasonal_effects_config,
                                         sim_factory.vector_world,
                                         sim_factory.clock, enable_ctypes)
    sim_factory.set_seasonal_effects_model(seasonal_effects)

    # Create health model
    health_model_class = config['health_model.__type__']
    health_model_config = config.subconfig('health_model')
    health_model = instantiate_class("pandemia.components.health_model",
                                     health_model_class, health_model_config, scale_factor,
                                     sim_factory.clock, enable_ctypes)
    sim_factory.set_health_model(health_model)

    # Create movement model
    movement_model_class = config['movement_model.__type__']
    movement_model_config = config.subconfig('movement_model')
    movement_model = instantiate_class("pandemia.components.movement_model",
                                       movement_model_class, movement_model_config, enable_ctypes)
    sim_factory.set_movement_model(movement_model)

    # Create hospitalization and death model
    hospitalization_and_death_model_class = config['hospitalization_and_death_model.__type__']
    hospitalization_and_death_model_config = config.subconfig('hospitalization_and_death_model')
    hospitalization_and_death_model =\
        instantiate_class("pandemia.components.hospitalization_and_death_model",
                          hospitalization_and_death_model_class,
                          hospitalization_and_death_model_config, enable_ctypes)
    sim_factory.set_hospitalization_and_death_model(hospitalization_and_death_model)

    # Create testing and contact tracing model
    testing_and_contact_tracing_model_class = config['testing_and_contact_tracing_model.__type__']
    testing_and_contact_tracing_model_config = config.subconfig('testing_and_contact_tracing_model')
    testing_and_contact_tracing_model =\
        instantiate_class("pandemia.components.testing_and_contact_tracing_model",
                          testing_and_contact_tracing_model_class,
                          testing_and_contact_tracing_model_config, enable_ctypes)
    sim_factory.set_testing_and_contact_tracing_model(testing_and_contact_tracing_model)

    # Create vaccination model
    vaccination_model_class = config['vaccination_model.__type__']
    vaccination_model_config = config.subconfig('vaccination_model')
    vaccination_model = instantiate_class("pandemia.components.vaccination_model",
                                          vaccination_model_class, vaccination_model_config,
                                          sim_factory.clock,
                                          health_model.number_of_strains,
                                          health_model.immunity_length,
                                          health_model.immunity_period_ticks,
                                          enable_ctypes)
    sim_factory.set_vaccination_model(vaccination_model)

    # Create regional mixing model
    regional_mixing_model_class = config['regional_mixing_model.__type__']
    regional_mixing_config = config.subconfig('regional_mixing_model')
    regional_mixing = instantiate_class("pandemia.components.regional_mixing_model",
                                        regional_mixing_model_class, regional_mixing_config,
                                        scale_factor,
                                        health_model.number_of_strains,
                                        sim_factory.vector_world,
                                        enable_ctypes)
    sim_factory.set_regional_mixing_model(regional_mixing)

    # Create input model
    input_class = config['input_model.__type__']
    input_config = config.subconfig('input_model')
    input = instantiate_class("pandemia.components.input_model",
                              input_class, input_config, scale_factor,
                              sim_factory.clock,
                              sim_factory.vector_world.number_of_regions,
                              vaccination_model.number_of_vaccines,
                              vaccination_model.age_groups,
                              enable_ctypes)
    sim_factory.set_input_model(input)

def get_historial_data():
    """Builds simulator and runs"""

    # Load simulation factory and build sim
    sim_factory = SimulationFactory.from_file(world_filepath)
    build_components(sim_factory)
    telemetry_bus = MessageBus()
    sim = sim_factory.new_sim(telemetry_bus)
    sim.input_model.new_input(None)
    sim.enable_parallel = False

    # Set validation parameter values
    sim.seasonal_effects_model.out_of_season_multiplier = 0.75
    sim.health_model.beta = 0.025
    sim.health_model.location_typ_multipliers['Square'] = 0.5
    sim.health_model.num_initial_infections_by_region['GB'] = [100000]
    sim.health_model.preexisting_sigma_multiplier = 0.5
    sim.health_model.preexisting_rho_multiplier = 0.5
    sim.regional_mixing_model.travel_multiplier = 10.0
    sim.input_model.max_transmission_control = 0.7
    sim.input_model.max_travel_control = 0.5

    sim.random_seed = 1

    # Run simulation and calculate error versus historical data
    sim.setup()
    sim.run()
    sim.calculate_error() # Scale historical deaths?

    return sim.average_deaths_dict_historical

def build_and_run(solution, solution_idx):
    """Builds simulator and runs"""

    # Load simulation factory and build sim
    sim_factory = SimulationFactory.from_file(world_filepath)
    build_components(sim_factory)
    telemetry_bus = MessageBus()
    sim = sim_factory.new_sim(telemetry_bus)
    sim.input_model.new_input(None)
    sim.enable_parallel = False

    # Set validation parameter values
    solution = np.divide(1, 1 + np.power(np.e, solution))

    sim.seasonal_effects_model.out_of_season_multiplier = (solution[0] * (1.00-0.25)) + 0.25
    sim.health_model.beta = (solution[1] * (0.05-0.01)) + 0.01
    sim.health_model.location_typ_multipliers['Square'] = (solution[2] * (1.0-0.0)) + 0.0
    sim.health_model.num_initial_infections_by_region['GB'] = [(solution[3] * (1000000-2000)) + 2000]
    sim.health_model.preexisting_sigma_multiplier = (solution[4] * (1.0-0.1)) + 0.1
    sim.health_model.preexisting_rho_multiplier = (solution[5] * (1.0-0.1)) + 0.1
    sim.regional_mixing_model.travel_multiplier = (solution[6] * (100.0-1.0)) + 1.0
    sim.input_model.max_transmission_control = (solution[7] * (1.0-0.25)) + 0.25
    sim.input_model.max_travel_control = (solution[8] * (1.0-0.25)) + 0.25

    sim.random_seed = 1

    # Run simulation and calculate error versus historical data
    sim.setup()
    sim.run()
    sim.average_deaths_dict_historical = copy.deepcopy(average_deaths_dict_historical)
    sim.calculate_error() # Scale historical deaths?

    return 1 / sim.error

def tabulate_results(num_configs):
    """Creates csv file of config ids matched to corresponding errors"""

    results_dict = {}
    for config_id in range(num_configs):
        with open('/tmp/validation_ga/configs/' + str(config_id) + '.csv', newline='') as csvfile:
            data = csv.reader(csvfile, delimiter=',')
            for row in data:
                if row[0] == 'result':
                    results_dict[config_id] = float(row[1])

    handle_results = open('/tmp/validation_ga/results.csv', 'w', newline='')
    writer_results = csv.writer(handle_results)
    for key in results_dict:
        row = [key, str(results_dict[key])]
        writer_results.writerow(row)
    handle_results.close()

world_filepath = 'fr.wld'

average_deaths_dict_historical = get_historial_data()

config_optimizer =\
{
    'num_generations': 20,
    'num_parents_mating': 5,
    'sol_per_pop': 10,
    'mutation_num_genes': 3,
    'num_genes': 9
}

fitness_func = build_and_run

def fitness_wrapper(solution):
    return fitness_func(solution, 0)

class PooledGA(pygad.GA):
    def cal_pop_fitness(self):
        global pool
        pop_fitness = pool.map(fitness_wrapper, self.population)
        pop_fitness = np.array(pop_fitness)
        return pop_fitness

def on_generation(ga_instance):
    """Prints generation number and cost of best solution"""
    print("Generation = {generation}".format(generation=ga_instance.generations_completed))
    print("Fitness    = {cost}".format(cost=\
          ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1]))

ga_instance = PooledGA(num_generations=config_optimizer['num_generations'],
                       num_parents_mating=config_optimizer['num_parents_mating'],
                       sol_per_pop=config_optimizer['sol_per_pop'],
                       mutation_num_genes=config_optimizer['mutation_num_genes'],
                       num_genes=config_optimizer['num_genes'],
                       fitness_func=fitness_func,
                       on_generation=on_generation,
                       gene_type=float)

if __name__ == "__main__":

    # Load the saved GA instance
    # filename = 'genetic'
    # ga_instance = pygad.load(filename=filename)

    # Run the genetic algorithm
    with Pool(processes=10) as pool:
        ga_instance.run()

    # Plot the fitness curve
    ga_instance.plot_fitness()

    # Return the details of the best solution
    solution, solution_fitness, solution_idx =\
        ga_instance.best_solution(ga_instance.last_generation_fitness)

    print("Best solution:")
    print(solution)

    # Save the GA instance
    filename = 'genetic'
    ga_instance.save(filename=filename)
