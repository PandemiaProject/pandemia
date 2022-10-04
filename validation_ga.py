
import pygad
import logging
import logging.config
from functools import partial
import csv
from multiprocessing import Pool
import numpy as np

from pandemia.messagebus import MessageBus
from pandemia.sim_factory import SimulationFactory

# Global module log
log = logging.getLogger()

def get_historial_data():
    """Builds simulator and runs"""

    # Load simulation factory and build sim
    sim_factory = SimulationFactory.from_file(world_filepath)
    sim_factory.build_components()
    telemetry_bus = MessageBus()
    sim = sim_factory.new_sim(telemetry_bus)
    sim.input_model.new_input(None)
    sim.enable_parallel = False

    # Set validation parameter values
    sim.seasonal_effects_model.out_of_season_multiplier = 0.75
    sim.health_model.beta = [0.025]
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

def build_and_run(solution, seed):
    """Builds simulator and runs"""

    # Load simulation factory and build sim
    sim_factory = SimulationFactory.from_file(world_filepath)
    sim_factory.build_components()
    telemetry_bus = MessageBus()
    sim = sim_factory.new_sim(telemetry_bus)
    sim.input_model.new_input(None)
    sim.enable_parallel = False

    # Set validation parameter values
    solution = np.divide(1, 1 + np.power(np.e, solution))

    sim.seasonal_effects_model.out_of_season_multiplier = (solution[0] * (1.00-0.5)) + 0.5
    sim.health_model.beta = [(solution[1] * (0.045-0.025)) + 0.025]
    sim.health_model.location_typ_multipliers['Square'] = (solution[2] * (0.5-0.1)) + 0.1
    sim.health_model.num_initial_infections_by_region['CN'] = [(solution[3] * (1000000-50000)) + 50000]
    sim.health_model.preexisting_sigma_multiplier = (solution[4] * (1.0-0.5)) + 0.5
    sim.health_model.preexisting_rho_multiplier = (solution[5] * (1.0-0.5)) + 0.5
    sim.regional_mixing_model.travel_multiplier = (solution[6] * (100.0-1.0)) + 1.0
    sim.input_model.max_transmission_control = (solution[7] * (1.0-0.5)) + 0.5
    sim.input_model.max_travel_control = (solution[8] * (1.0-0.5)) + 0.5
    sim.regional_mixing_model.interpolation = (solution[9] * (0.01-0.0)) + 0.0

    sim.random_seed = seed

    # Run simulation and calculate error versus historical data
    sim.setup()
    sim.run()
    # sim.average_deaths_dict_historical = copy.deepcopy(average_deaths_dict_historical)
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

seed = 0

world_filepath = 'Scenarios\Global_Grid\global_grid_world.wld'

# average_deaths_dict_historical = get_historial_data()

config_optimizer =\
{
    'num_generations': 100,
    'num_parents_mating': 4,
    'sol_per_pop': 8,
    'mutation_num_genes': 2,
    'num_genes': 10
}

def fitness_func(solution, solution_idx):
    return build_and_run(solution, 0)

class PooledGA(pygad.GA):
    def cal_pop_fitness(self):
        global pool
        seed = self.generations_completed
        pop_fitness = pool.map(partial(build_and_run, seed=seed), self.population)
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
    with Pool(processes=8) as pool:
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
