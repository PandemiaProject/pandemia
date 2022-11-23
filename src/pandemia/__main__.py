"""Pandemia is an individual-based pandemic simulator."""

import os.path as osp
import sys
import logging
import logging.config
import time
import cProfile, pstats

from .messagebus import MessageBus
from .sim_factory import SimulationFactory

from .version import VERSION
from .config import Config

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

    policy = None
    sim.input_model.new_input(policy)

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

    sim.calculate_cost(policy)

if __name__ == "__main__":
    main()