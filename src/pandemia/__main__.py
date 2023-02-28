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

    ################################################################################################

    # # sim_factory.clock.simulation_length_days = 10 # 730 # can manually adjust the sim length here
    # random_seed = 5
    # enable_parallel = False
    # num_jobs = 1

    # # sim_factory.config['health_model']['num_initial_infections_by_region_by_strain'] = {'SA': [500000]}

    # ### 1.0 Transmissibility ###
    # # sim_factory.config['health_model']['beta'] = [0.035]
    # ### 2.0 Transmissibility ###
    # sim_factory.config['health_model']['beta'] = [0.070]

    # ### 0.5% IFR ###
    # sim_factory.config['health_model']['preset_weights_by_age'] = {
    #     0:  {'preset_0': 1.00000, 'preset_1': 0.00000},
    #     5:  {'preset_0': 0.99988, 'preset_1': 0.00012},
    #     10: {'preset_0': 0.99988, 'preset_1': 0.00012},
    #     15: {'preset_0': 0.99976, 'preset_1': 0.00024},
    #     20: {'preset_0': 0.99976, 'preset_1': 0.00024},
    #     25: {'preset_0': 0.99952, 'preset_1': 0.00048},
    #     30: {'preset_0': 0.99928, 'preset_1': 0.00072},
    #     35: {'preset_0': 0.99892, 'preset_1': 0.00108},
    #     40: {'preset_0': 0.99820, 'preset_1': 0.00180},
    #     45: {'preset_0': 0.99724, 'preset_1': 0.00276},
    #     50: {'preset_0': 0.99568, 'preset_1': 0.00432},
    #     55: {'preset_0': 0.99316, 'preset_1': 0.00684},
    #     60: {'preset_0': 0.98932, 'preset_1': 0.01068},
    #     65: {'preset_0': 0.98332, 'preset_1': 0.01668},
    #     70: {'preset_0': 0.97396, 'preset_1': 0.02604},
    #     75: {'preset_0': 0.95932, 'preset_1': 0.04068},
    #     80: {'preset_0': 0.93640, 'preset_1': 0.06360},
    #     85: {'preset_0': 0.90064, 'preset_1': 0.09936},
    #     90: {'preset_0': 0.80572, 'preset_1': 0.19428}
    # }
    # ### 1.5% IFR ###
    # # sim_factory.config['health_model']['preset_weights_by_age'] = {
    # #     0:  {'preset_0': 1.00000, 'preset_1': 0.00000},
    # #     5:  {'preset_0': 0.99964, 'preset_1': 0.00036},
    # #     10: {'preset_0': 0.99964, 'preset_1': 0.00036},
    # #     15: {'preset_0': 0.99928, 'preset_1': 0.00072},
    # #     20: {'preset_0': 0.99928, 'preset_1': 0.00072},
    # #     25: {'preset_0': 0.99856, 'preset_1': 0.00144},
    # #     30: {'preset_0': 0.99784, 'preset_1': 0.00216},
    # #     35: {'preset_0': 0.99676, 'preset_1': 0.00324},
    # #     40: {'preset_0': 0.99460, 'preset_1': 0.00540},
    # #     45: {'preset_0': 0.99172, 'preset_1': 0.00828},
    # #     50: {'preset_0': 0.98704, 'preset_1': 0.01296},
    # #     55: {'preset_0': 0.97948, 'preset_1': 0.02052},
    # #     60: {'preset_0': 0.96796, 'preset_1': 0.03204},
    # #     65: {'preset_0': 0.94996, 'preset_1': 0.05004},
    # #     70: {'preset_0': 0.92188, 'preset_1': 0.07812},
    # #     75: {'preset_0': 0.87796, 'preset_1': 0.12204},
    # #     80: {'preset_0': 0.80920, 'preset_1': 0.19080},
    # #     85: {'preset_0': 0.70192, 'preset_1': 0.29808},
    # #     90: {'preset_0': 0.41716, 'preset_1': 0.58284}
    # # }

    # sim_factory.config['reporters']['csv.StrainCounts']['filename'] = "/tmp/strain_counts_heterogeneous_" + str(random_seed) + ".csv"
    # sim_factory.config['reporters']['csv.DeathCounts']['filename'] = "/tmp/death_counts_heterogeneous_" + str(random_seed) + ".csv"
    # sim_factory.config['reporters']['plot.PlotInfected']['filename'] = "/tmp/infected_heterogeneous_" + str(random_seed) + ".png"
    # sim_factory.config['reporters']['plot.PlotDeaths']['daily_deaths_filename'] = "/tmp/daily_deaths_heterogeneous_" + str(random_seed) + ".png"
    # sim_factory.config['reporters']['plot.PlotDeaths']['cumulative_deaths_filename'] = "/tmp/cumulative_deaths_heterogeneous_" + str(random_seed) + ".png"
    # sim_factory.config['reporters']['plot.PlotDeaths']['deaths_by_country_filename'] = "/tmp/deaths_by_country_" + str(random_seed) + ".csv"

    ################################################################################################

    # Build simulation components
    sim_factory.build_components()

    # Build reporters according to config
    telemetry_bus = MessageBus()
    sim_factory.build_reporters(telemetry_bus)

    sim = sim_factory.new_sim(telemetry_bus)

    ################################################################################################

    # sim.random_seed = random_seed
    # sim.enable_parallel = enable_parallel
    # sim.num_jobs= num_jobs

    ################################################################################################

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