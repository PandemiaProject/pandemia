"""Model validation tool based on grid search (work in progress)"""

from copy import deepcopy
import logging
import logging.config
import csv

from joblib import Parallel, delayed

from .messagebus import MessageBus
from .sim_factory import SimulationFactory

# Global module log
log = logging.getLogger()

def build_and_run(validation_config, num_configs, world_filepath):
    """Builds simulator and runs"""

    config_id = validation_configs.index(validation_config)

    # Load simulation factory and build sim
    sim_factory = SimulationFactory.from_file(world_filepath)
    sim_factory.build_components()
    telemetry_bus = MessageBus()
    sim = sim_factory.new_sim(telemetry_bus)
    sim.input_model.new_input(None)
    sim.enable_parallel = False

    # set validation parameter values
    sim.random_seed =\
        validation_config['random_seed']
    sim.seasonal_effects_model.out_of_season_multiplier =\
        validation_config['out_of_season_multiplier']
    sim.health_model.beta =\
        [validation_config['beta']]
    sim.health_model.location_typ_multipliers['Square'] =\
        validation_config['location_typ_multipliers[Square]']
    sim.health_model.num_initial_infections_by_region['GB'] =\
        [validation_config['num_initial_infections']]
    sim.health_model.preexisting_sigma_multiplier =\
        validation_config['preexisting_sigma_multiplier']
    sim.health_model.preexisting_rho_multiplier =\
        validation_config['preexisting_rho_multiplier']
    sim.travel_model.travel_multiplier =\
        validation_config['travel_multiplier']
    sim.input_model.max_transmission_control =\
        validation_config['max_transmission_control']
    sim.input_model.max_travel_control =\
        validation_config['max_travel_control']

    # Run simulation and calculate error versus historical data
    sim.setup()
    sim.run()
    sim.calculate_error() # Scale historical deaths?

    # Save results of simulation
    handle_config = open('/tmp/validation_gs/configs/' + str(config_id) + '.csv', 'w', newline='')
    writer_config = csv.writer(handle_config)
    for key in validation_config:
        row = [key, str(validation_config[key])]
        writer_config.writerow(row)
    writer_config.writerow(['result', str(sim.error)])
    handle_config.close()

    print("Completed simulation for config: ", config_id, "/", num_configs)

def tabulate_results(num_configs):
    """Creates csv file of config ids matched to corresponding errors"""

    results_dict = {}
    for config_id in range(num_configs):
        with open('/tmp/validation_gs/configs/' + str(config_id) + '.csv', newline='') as csvfile:
            data = csv.reader(csvfile, delimiter=',')
            for row in data:
                if row[0] == 'result':
                    results_dict[config_id] = float(row[1])

    handle_results = open('/tmp/validation_gs/results.csv', 'w', newline='')
    writer_results = csv.writer(handle_results)
    for key in results_dict:
        row = [key, str(results_dict[key])]
        writer_results.writerow(row)
    handle_results.close()

validation_config_default =\
{
    'random_seed':                      1,
    'out_of_season_multiplier':         0.75,
    'beta':                             0.025,
    'location_typ_multipliers[Square]': 0.5,
    'num_initial_infections':           100000,
    'preexisting_sigma_multiplier':     0.5,
    'preexisting_rho_multiplier':       0.5,
    'travel_multiplier':                10.0,
    'max_transmission_control':         0.7,
    'max_travel_control':               0.5
}

def generate_validation_configs(validation_config_default):
    """Generate sample configs"""

    validation_configs = []
    for b in range(5):
        for o in range(3):
            new_config = deepcopy(validation_config_default)
            new_config['beta'] = 0.02 + (b * 0.001)
            new_config['out_of_season_multiplier'] = 0.7 + (o * 0.01)
            validation_configs.append(new_config)

    return validation_configs

validation_configs = generate_validation_configs(validation_config_default)
num_configs = len(validation_configs)

world_filepath = 'gb.wld'

Parallel(n_jobs=15)(delayed(build_and_run)(validation_config, num_configs, world_filepath)
                for validation_config in validation_configs)

tabulate_results(num_configs)
