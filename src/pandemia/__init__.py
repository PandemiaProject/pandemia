"""Pandemia is an individual-based multi-region epidemic simulator."""

import os.path as osp
import sys
import logging
import logging.config
import time
import cProfile, pstats

from pandemia.messagebus import MessageBus
from pandemia.sim_factory import SimulationFactory

from pandemia.version import VERSION
from pandemia.config import Config

# Global module log
log = logging.getLogger()


def main():
    """Main pandemia entry point. Two command-line arguments may be passed, the first is required
    while the second is optional. The first should be the filepath to a simulation config, while the
    second should be a filepath for either loading or saving a simulation factory. Building a
    simulation world can take a while, hence the loading and saving feature."""
    print(f"pandemia {VERSION}")

    # A configuration file must be given before building a simulation. This specifices all the
    # parameters for the simulation.
    try:
        assert len(sys.argv) > 1
    except AssertionError:
        log.warning("No config given")
        exit(1)

    # System configuration and setup
    if len(sys.argv) > 2 and osp.isfile(sys.argv[2]):
        sim_factory = SimulationFactory.from_file(sys.argv[2])
        logging.config.dictConfig(sim_factory.config['logging'])
        log.warning("Existing factory loaded from %s", sys.argv[2])
    else:
        config = Config(sys.argv[1])
        sim_factory = SimulationFactory(config)
        logging.config.dictConfig(sim_factory.config['logging'])

        log.info("State info:")
        log.info("  Run ID: %s", sim_factory.run_id)
        log.info("  pandemia version: %s", sim_factory.pandemia_version)
        log.info("  Created at: %s", sim_factory.created_at)

        sim_factory.build_clock_and_world()

        if len(sys.argv) > 2:
            log.info("Writing to state file: %s", sys.argv[2])
            sim_factory.to_file(sys.argv[2])

    # sim_factory.clock.simulation_length_days = 730 # can manually adjust the sim length here

    # Build simulation components
    sim_factory.build_components()

    # Build reporters according to config
    telemetry_bus = MessageBus()
    sim_factory.build_reporters(telemetry_bus)

    sim = sim_factory.new_sim(telemetry_bus)

    # ############## Input ##############

    input_arrays = None
    sim.input_model.new_input(input_arrays)

    # ############## Manual Parameter Adjustment ##############

    # import numpy as np
    # sim.enable_parallel = False
    # sim.random_seed = 1
    # sim.seasonal_effects_model.out_of_season_multiplier = 0.75
    # sim.health_model.beta = [0.035]
    # mutation_matrix = [[1.0]]
    # sim.health_model.mutation_matrix = np.asarray(mutation_matrix, dtype=float).flatten()
    # sim.health_model.location_typ_multipliers['Square'] = 0.2
    # sim.health_model.num_initial_infections_by_region['CN'] = [500000]
    # sim.health_model.preexisting_sigma_multiplier = 0.5
    # sim.health_model.preexisting_rho_multiplier = 0.5
    # sim.regional_mixing_model.travel_multiplier = 600
    # sim.input_model.max_transmission_control = 0.9
    # sim.input_model.max_travel_control = 0.9
    # sim.regional_mixing_model.interpolation = 0.00

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
