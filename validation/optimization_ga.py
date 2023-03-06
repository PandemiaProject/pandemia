"""Policy optimization tool based on genetic algorithms (work in progress)"""

import pygad
import logging
import logging.config
from functools import partial
from multiprocessing import Pool
import numpy as np

from .messagebus import MessageBus
from .sim_factory import SimulationFactory

# Global module log
log = logging.getLogger()

class Policy():
    """Represents a policy"""

    def __init__(self, T, R, A, V):

        self.lockdown_policy            = np.full((T, R), -1, dtype=np.int32)
        self.border_closure_policy      = np.full((T, R), 1.0, dtype=np.float64)
        self.facemask_policy            = np.full((T, R), -1, dtype=np.int32)
        self.random_testing_policy      = np.full((T, R), 0, dtype=np.int32)
        self.symptomatic_testing_policy = np.full((T, R), 0, dtype=np.int32)
        self.contact_testing_policy     = np.full((T, R), 0, dtype=np.int32)
        self.vaccination_policy         = np.full((T, R, A, V), 0, dtype=np.int32)

def fitness_func(solution, solution_idx):
    """Builds simulator and runs"""

    # Create new policy

    policy = Policy(num_days, num_regions, num_age_groups, num_vaccines)

    solution = np.divide(1, 1 + np.power(np.e, solution))

    TOTAL_DOSES = 10000000 # The total number of doses available
    VACCINATION_RATE = 0.006 # The maximum proportion of population a country can vaccinate each day
    LEN_SUPPLY_PERIOD = 30   # The time period is divided into supply periods of this length in days

    vaccination_rate = np.full((num_regions), VACCINATION_RATE, dtype=np.float64)
    num_supply_periods = (num_days // LEN_SUPPLY_PERIOD) + 1
    supply = np.full((num_supply_periods), int(TOTAL_DOSES / num_supply_periods), dtype=np.int32)

    solution_part_1 = solution[0: num_regions * num_supply_periods]
    solution_part_2 = solution[num_regions * num_supply_periods:]

    vac_sol = np.absolute(solution_part_1.reshape((num_supply_periods, num_regions)))
    vac_sol = ((vac_sol.T / np.sum(vac_sol, axis=1)) * supply).T

    dist = np.absolute(solution_part_2)
    dist = np.reshape(dist, (num_supply_periods, num_regions, num_age_groups))
    for l in range(num_supply_periods):
        dist_sum = dist[l].sum(axis=1)
        dist_sum[dist_sum == 0] = 1
        dist[l] = np.transpose(np.transpose(dist[l]) / dist_sum)

    num_can_vaccinate_each_day = (np.multiply(vaccination_rate, population_sizes)).astype(int)

    vaccination_sol = np.zeros((num_days, num_regions, num_age_groups), dtype=np.int32)
    for day in range(num_days):
        supply_period = day // LEN_SUPPLY_PERIOD
        share = np.minimum((vac_sol[supply_period] / LEN_SUPPLY_PERIOD).astype(int),
                            num_can_vaccinate_each_day)
        vaccination_sol[day] = (np.multiply(share[:, None], dist[supply_period])).astype(int)

    policy.vaccination_policy = vaccination_sol[:, :, :, np.newaxis]

    average_cost = 0
    num_valid_runs = 0
    costs = []

    for seed in seeds:

        # Build sim
        sim_factory = SimulationFactory.from_file(world_filepath)
        sim_factory.build_components()
        telemetry_bus = MessageBus()
        sim = sim_factory.new_sim(telemetry_bus)
        sim.enable_parallel = False
        sim.random_seed = seed

        # Run sim
        sim.policy_maker_model.new_policy(policy)
        sim.setup()
        sim.run()
        cost = sim.calculate_cost(policy)

        # Calculate cost
        if cost != 0:
            average_cost += cost
            costs.append(cost)
            num_valid_runs += 1
    
    if num_valid_runs != 0:
        average_cost /= num_valid_runs
    else:
        average_cost = 9999999999

    print(costs, average_cost)

    return 1 / average_cost

seed = 0

world_filepath = 'gb_world.wld'

sim_factory = SimulationFactory.from_file(world_filepath)
sim_factory.build_components()
telemetry_bus = MessageBus()
sim = sim_factory.new_sim(telemetry_bus)

num_days = sim.policy_maker_model.simulation_length_days
num_regions = sim.policy_maker_model.number_of_regions
num_vaccines = sim.policy_maker_model.number_of_vaccines
num_age_groups = sim.policy_maker_model.number_of_vaccination_age_groups

assert num_vaccines == 1

scale_factor = sim.vector_world.scale_factor
rescale_factor = 1 / scale_factor
population_sizes =\
    np.array([vr.number_of_agents * rescale_factor for vr in sim.vector_regions], dtype=np.float64)

seeds = [0, 1, 2, 3, 4]

config_optimizer =\
{
    'num_generations': 100,
    'num_parents_mating': 6,
    'sol_per_pop': 12,
    'num_genes': 13 + 39
}

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
                       num_genes=config_optimizer['num_genes'],
                       fitness_func=fitness_func,
                       on_generation=on_generation,
                       gene_type=float)

if __name__ == "__main__":

    # Load the saved GA instance
    # filename = 'genetic'
    # ga_instance = pygad.load(filename=filename)

    # Run the genetic algorithm
    with Pool(processes=12) as pool:
        ga_instance.run()

    # Plot the fitness curve
    ga_instance.plot_fitness()

    # Return the details of the best solution
    solution, solution_fitness, solution_idx =\
        ga_instance.best_solution(ga_instance.last_generation_fitness)

    # Save the GA instance
    filename = 'genetic'
    ga_instance.save(filename=filename)
