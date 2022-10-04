"""Pandemia is an individual-based multi-region epidemic simulator."""

import os.path as osp
import sys
import logging
import logging.config
import time
import cProfile, pstats

from pandemia.utils import instantiate_class
from pandemia.messagebus import MessageBus
from pandemia.clock import Clock
from pandemia.sim_factory import SimulationFactory

from pandemia.version import VERSION
from pandemia.config import Config

# Global module log
log = logging.getLogger()

def build_clock_and_world(sim_factory):
    """Builds the model on which the simulator acts. The model consists of a clock, to represent
    time, and a world, representing regions. Each region consists of agents and locations, with
    agents able to perform activities in selected locations. There is also a scale factor, which can
    be used to rescale the numbers of agents and locations in each region"""

    config = sim_factory.config

    scale_factor = config['scale_factor']

    # Create clock
    _clock = Clock(config['tick_length_s'], config['simulation_length_days'], config['epoch'])
    log.info("New clock created at %s, tick_length = %i, simulation_days = %i, week_offset = %i",
             _clock.epoch, _clock.tick_length_s, _clock.simulation_length_days,
             _clock.epoch_week_offset)
    sim_factory.set_clock(_clock)

    # Create world
    world_factory_class = config['world_factory.__type__']
    world_factory_config = config.subconfig('world_factory')
    world_factory = instantiate_class("pandemia.world.world_factory", world_factory_class,
                                      world_factory_config, sim_factory.clock, scale_factor)
    world = world_factory.get_world()
    vector_world = world.vectorize_world()
    sim_factory.set_vector_world(vector_world)

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

def build_reporters(telemetry_bus, config):
    """Instantiates reporters, which record output data on the simulation for analysis"""

    if config['reporters'] is not None:
        for reporter_class, reporter_config in config['reporters'].items():
            log.info(f"Creating reporter '{reporter_class}'...")
            instantiate_class("pandemia.reporters", reporter_class, telemetry_bus,
                              Config(_dict=reporter_config))

def main():
    """Main pandemia entry point"""
    print(f"pandemia {VERSION}")

    try:
        assert len(sys.argv) > 1
    except AssertionError:
        log.warning("No config given")
        exit(1)

    # System config/setup
    if len(sys.argv) > 2 and osp.isfile(sys.argv[2]):
        sim_factory = SimulationFactory.from_file(sys.argv[2])
        logging.config.dictConfig(sim_factory.config['logging'])
        log.warning("Existing factory loaded from %s", sys.argv[2])
    else:
        config = Config(sys.argv[1])
        sim_factory = SimulationFactory(config)
        logging.config.dictConfig(sim_factory.config['logging'])

        # Summarise the sim_factory
        log.info("State info:")
        log.info("  Run ID: %s", sim_factory.run_id)
        log.info("  pandemia version: %s", sim_factory.pandemia_version)
        log.info("  Created at: %s", sim_factory.created_at)

        build_clock_and_world(sim_factory)

        # If a second parameter is given, use this for the statefile name
        if len(sys.argv) > 2:
            log.info("Writing to state file: %s", sys.argv[2])
            sim_factory.to_file(sys.argv[2])

    # sim_factory.clock.simulation_length_days = 730 # 365

    # Build simulation components
    build_components(sim_factory)

    # Build reporters according to config
    telemetry_bus = MessageBus()
    build_reporters(telemetry_bus, sim_factory.config)

    sim = sim_factory.new_sim(telemetry_bus)

    # ############## Input ##############

    input_arrays = None
    sim.input_model.new_input(input_arrays)

    # ############## set validation parameter values (move all this etc...) ##############
    # sim.enable_parallel = False
    # sim.random_seed = 13
    # import numpy as np
    # solution = np.array([-3.6723102, -0.64972947, 2.72811164, 4.44037522, 3.13302174, -2.28823345,
    #                      -2.06260188, 2.02718859, 2.99305056, 3.52555033])

    # solution = np.divide(1, 1 + np.power(np.e, solution))

    # sim.seasonal_effects_model.out_of_season_multiplier = (solution[0] * (1.00-0.25)) + 0.25
    # sim.health_model.beta = [(solution[1] * (0.05-0.01)) + 0.01]
    # sim.health_model.location_typ_multipliers['Square'] = (solution[2] * (1.0-0.0)) + 0.0
    # sim.health_model.num_initial_infections_by_region['CN'] = [(solution[3] * (1000000-2000)) + 2000]
    # sim.health_model.preexisting_sigma_multiplier = (solution[4] * (1.0-0.1)) + 0.1
    # sim.health_model.preexisting_rho_multiplier = (solution[5] * (1.0-0.1)) + 0.1
    # sim.regional_mixing_model.travel_multiplier = (solution[6] * (100.0-1.0)) + 1.0
    # sim.input_model.max_transmission_control = (solution[7] * (1.0-0.25)) + 0.25
    # sim.input_model.max_travel_control = (solution[8] * (1.0-0.25)) + 0.25
    # sim.regional_mixing_model.interpolation = (solution[9] * (0.25-0.0)) + 0.0

    # sim.seasonal_effects_model.out_of_season_multiplier = 0.75
    # sim.health_model.beta = [0.035]
    # import numpy as np
    # mutation_matrix = [[0.9995, 0.0005], [0.0, 1.0]]
    # sim.health_model.mutation_matrix = np.asarray(mutation_matrix, dtype=float).flatten()
    # sim.health_model.location_typ_multipliers['Square'] = 0.2
    # sim.health_model.num_initial_infections_by_region['CN'] = [500000]
    # sim.health_model.preexisting_sigma_multiplier = 0.5
    # sim.health_model.preexisting_rho_multiplier = 0.5
    # sim.regional_mixing_model.travel_multiplier = 600
    # sim.input_model.max_transmission_control = 0.9
    # sim.input_model.max_travel_control = 0.9
    # sim.regional_mixing_model.interpolation = 0.00

    # print("out_of_season_multiplier: ", sim.seasonal_effects_model.out_of_season_multiplier)
    # print("beta: ", sim.health_model.beta)
    # print("location_typ_multipliers['Square']: ", sim.health_model.location_typ_multipliers['Square'])
    # print("num_initial_infections_by_region['CN']: ", sim.health_model.num_initial_infections_by_region['CN'][0])
    # print("preexisting_sigma_multiplier: ", sim.health_model.preexisting_sigma_multiplier)
    # print("preexisting_rho_multiplier: ", sim.health_model.preexisting_rho_multiplier)
    # print("travel_multiplier: ", sim.regional_mixing_model.travel_multiplier)
    # print("max_transmission_control: ", sim.input_model.max_transmission_control)
    # print("max_travel_control: ", sim.input_model.max_travel_control)
    # print("interpolation: ", sim.regional_mixing_model.interpolation)

    # ############## Setup ##############

    sim.setup()

    # ############### Run ###############

    if sim_factory.config['profile']:
        profiler = cProfile.Profile()
        profiler.enable()

    t0 = time.time()

    sim.run()

    t1 = time.time()
    total_time = t1 - t0

    if sim_factory.config['profile']:
        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats('tottime')
        stats.strip_dirs()
        stats.dump_stats('profiler_output')
        file = open(str(sim_factory.run_id) + '.csv', 'w')
        profile = pstats.Stats('.\profiler_output', stream=file)
        profile.sort_stats('tottime')
        profile.print_stats(sim_factory.config['pstats_records'])
        file.close()

    log.info("Simulation finished successfully in " + str(round(total_time, 2)) + " seconds.")

    # ############## Output ##############

    sim.calculate_cost(input_arrays)
