"""Pandemia is an individual-based pandemic simulator."""

import os.path as osp
import sys
import logging
import logging.config
import time
import cProfile, pstats
from .utils import instantiate_class

from .messagebus import MessageBus
from .sim_factory import SimulationFactory
from .world import VectorWorld

from .version import VERSION
from .config import Config

# Global module log
log = logging.getLogger()

def build_world(config):
    """Builds the world on which the simulator acts. The world consists of regions. Each region
    consists of agents and locations, with agents able to perform activities in selected locations.
    There is also a scale factor, which can be used to rescale the numbers of agents and locations
    in each region"""

    scale_factor = config['scale_factor']

    # Create world
    world_factory_class = config['world_factory.__type__']
    world_factory_config = config.subconfig('world_factory')
    world_factory = instantiate_class("pandemia.world.world_factory", world_factory_class,
                                      world_factory_config, scale_factor)
    world = world_factory.get_world()
    vector_world = world.vectorize_world()

    return vector_world

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

    # Load or build world
    config = Config(sys.argv[1])
    logging.config.dictConfig(config['logging'])
    if len(sys.argv) > 2 and osp.isfile(sys.argv[2]):
        log.info('Reading data from %s...', sys.argv[2])
        vector_world = VectorWorld.from_file(sys.argv[2])
        log.warning("Existing world loaded from %s", sys.argv[2])
    else:
        vector_world = build_world(config)
        if len(sys.argv) > 2:
            log.info("Writing to world file: %s", sys.argv[2])
            vector_world.to_file(sys.argv[2])

    # Instantiate simulation factory
    sim_factory = SimulationFactory(config, vector_world)

    # Build simulation clock
    sim_factory.build_clock()

    # Build simulation components
    sim_factory.build_components()

    # Build reporters according to config
    telemetry_bus = MessageBus()
    sim_factory.build_reporters(telemetry_bus)

    sim = sim_factory.new_sim(telemetry_bus)

    # ############## Policy Maker ##############

    policy = None
    sim.policy_maker_model.new_policy(policy)

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