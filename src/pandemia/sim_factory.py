"""Builds a Simulator object. The SimulationFactory object builds a Simulator object gradually, by
adding the various model layers and components."""

# Allows classes to return their own type, e.g. from_file below
from __future__ import annotations

import logging
import uuid
import pickle
from datetime import datetime
from .utils import instantiate_class

from .config import Config
from .clock import Clock
from .version import VERSION
from .world import VectorWorld
from .simulator import Simulator
from .components.seasonal_effects_model import SeasonalEffectsModel
from .components.health_model import HealthModel
from .components.movement_model import MovementModel
from .components.hospitalization_and_death_model import HospitalizationAndDeathModel
from .components.testing_and_contact_tracing_model import TestingAndContactTracingModel
from .components.vaccination_model import VaccinationModel
from .components.travel_model import TravelModel
from .components.input_model import InputModel

log = logging.getLogger("sim_state")

class SimulationFactory:
    """Class that allows for gradual composition of a number of components, eventually outputting
    a Simulator object that can be used to run simulations with the given config."""

    def __init__(self, config: Config):

        # Static info
        self.pandemia_version = VERSION
        self.created_at       = datetime.now()
        self.run_id           = uuid.uuid4().hex

        log.info("Simulation state created at %s with ID=%s", self.created_at, self.run_id)

        # Config
        self.config = config

        # Clock
        self.clock = None

        # World
        self.vector_world = None

        # Components
        self.input_model = None
        self.seasonal_effects_model = None
        self.health_model = None
        self.movement_model = None
        self.hospitalization_and_death_model = None
        self.testing_and_contact_tracing_model = None
        self.vaccination_model = None
        self.travel_model = None

    def set_clock(self, clock: Clock) -> None:
        """Sets clock"""
        self.clock = clock

    def set_vector_world(self, vector_world: VectorWorld) -> None:
        """Sets world model"""
        self.vector_world = vector_world

    def set_seasonal_effects_model(self, seasonal_effects_model: SeasonalEffectsModel) -> None:
        """Sets seasonal effects model"""
        self.seasonal_effects_model = seasonal_effects_model

    def set_health_model(self, health_model: HealthModel) -> None:
        """Sets health model"""
        self.health_model = health_model

    def set_movement_model(self, movement_model: MovementModel) -> None:
        """Sets movement model"""
        self.movement_model = movement_model

    def set_hospitalization_and_death_model(self, hospitalization_and_death_model:
                                            HospitalizationAndDeathModel) -> None:
        """Sets hospitalization and death model"""
        self.hospitalization_and_death_model = hospitalization_and_death_model

    def set_testing_and_contact_tracing_model(self, testing_and_contact_tracing_model:
                                              TestingAndContactTracingModel) -> None:
        """Sets testing and contact tracing model"""
        self.testing_and_contact_tracing_model = testing_and_contact_tracing_model

    def set_vaccination_model(self, vaccination_model: VaccinationModel) -> None:
        """Sets vaccination model"""
        self.vaccination_model = vaccination_model

    def set_travel_model(self, travel_model: TravelModel) -> None:
        """Sets regional mixing model"""
        self.travel_model = travel_model

    def set_input_model(self, input_model: InputModel) -> None:
        """Sets seasonal effects model"""
        self.input_model = input_model

    def build_clock_and_world(self):
        """Builds the model on which the simulator acts. The model consists of a clock, to represent
        time, and a world, representing regions. Each region consists of agents and locations, with
        agents able to perform activities in selected locations. There is also a scale factor, which
        can be used to rescale the numbers of agents and locations in each region"""

        config = self.config

        scale_factor = config['scale_factor']

        # Create clock
        _clock = Clock(config['tick_length_s'], config['simulation_length_days'], config['epoch'])
        log.info("New clock created at %s, tick_length = %i, simulation_days = %i",
                 _clock.epoch, _clock.tick_length_s, _clock.simulation_length_days)
        self.set_clock(_clock)

        # Create world
        world_factory_class = config['world_factory.__type__']
        world_factory_config = config.subconfig('world_factory')
        world_factory = instantiate_class("pandemia.world.world_factory", world_factory_class,
                                        world_factory_config, self.clock, scale_factor)
        world = world_factory.get_world()
        vector_world = world.vectorize_world()
        self.set_vector_world(vector_world)

    def build_components(self):
        """The model also features a number of components, representing agent mobility, health and
        number of other features. There is also a scale factor, which can be used to rescale the
        numbers of tests and vaccine doses given as input, the numbers of agents travelling between
        the regions, and the final output."""

        config = self.config

        scale_factor = config['scale_factor']

        # Create seasonal effects model
        seasonal_effects_class = config['seasonal_effects_model.__type__']
        seasonal_effects_config = config.subconfig('seasonal_effects_model')
        seasonal_effects = instantiate_class("pandemia.components.seasonal_effects_model",
                                            seasonal_effects_class, seasonal_effects_config,
                                            self.vector_world,
                                            self.clock)
        self.set_seasonal_effects_model(seasonal_effects)

        # Create health model
        health_model_class = config['health_model.__type__']
        health_model_config = config.subconfig('health_model')
        health_model = instantiate_class("health_model",
                                        health_model_class, health_model_config, scale_factor,
                                        self.clock)
        self.set_health_model(health_model)

        # Create movement model
        movement_model_class = config['movement_model.__type__']
        movement_model_config = config.subconfig('movement_model')
        movement_model = instantiate_class("movement_model",
                                        movement_model_class, movement_model_config)
        self.set_movement_model(movement_model)

        # Create hospitalization and death model
        hospitalization_and_death_model_class = config['hospitalization_and_death_model.__type__']
        hospitalization_and_death_model_config = config.subconfig('hospitalization_and_death_model')
        hospitalization_and_death_model =\
            instantiate_class("hospitalization_and_death_model",
                            hospitalization_and_death_model_class,
                            hospitalization_and_death_model_config)
        self.set_hospitalization_and_death_model(hospitalization_and_death_model)

        # Create testing and contact tracing model
        testing_and_contact_tracing_model_class =\
            config['testing_and_contact_tracing_model.__type__']
        testing_and_contact_tracing_model_config =\
            config.subconfig('testing_and_contact_tracing_model')
        testing_and_contact_tracing_model =\
            instantiate_class(".testing_and_contact_tracing_model",
                            testing_and_contact_tracing_model_class,
                            testing_and_contact_tracing_model_config)
        self.set_testing_and_contact_tracing_model(testing_and_contact_tracing_model)

        # Create vaccination model
        vaccination_model_class = config['vaccination_model.__type__']
        vaccination_model_config = config.subconfig('vaccination_model')
        vaccination_model = instantiate_class(".vaccination_model",
                                            vaccination_model_class, vaccination_model_config,
                                            self.clock,
                                            health_model.number_of_strains,
                                            health_model.immunity_length,
                                            health_model.immunity_period_ticks)
        self.set_vaccination_model(vaccination_model)

        # Create regional mixing model
        travel_model_class = config['travel_model.__type__']
        travel_config = config.subconfig('travel_model')
        travel = instantiate_class(".travel_model",
                                            travel_model_class, travel_config,
                                            scale_factor,
                                            health_model.number_of_strains,
                                            self.vector_world)
        self.set_travel_model(travel)

        # Create input model
        input_class = config['input_model.__type__']
        input_config = config.subconfig('input_model')
        input = instantiate_class(".input_model",
                                input_class, input_config, scale_factor,
                                self.clock,
                                self.vector_world.number_of_regions,
                                vaccination_model.number_of_vaccines,
                                vaccination_model.age_groups)
        self.set_input_model(input)

    def build_reporters(self, telemetry_bus):
        """Instantiates reporters, which record output data on the simulation for analysis"""

        config = self.config

        if config['reporters'] is not None:
            for reporter_class, reporter_config in config['reporters'].items():
                log.info(f"Creating reporter '{reporter_class}'...")
                instantiate_class("pandemia.reporters", reporter_class, telemetry_bus,
                                Config(_dict=reporter_config))

    def new_sim(self, telemetry_bus):
        """Return a new simulator based on the given config.

        Telemetry data will be sent to the telemetry_bus provided (of type MessageBus)
        """

        if self.vector_world is None:
            raise ValueError("No world defined.")

        sim = Simulator(self.config,
                        self.clock,
                        self.vector_world,
                        self.seasonal_effects_model,
                        self.health_model,
                        self.movement_model,
                        self.hospitalization_and_death_model,
                        self.testing_and_contact_tracing_model,
                        self.vaccination_model,
                        self.travel_model,
                        self.input_model,
                        telemetry_bus)

        return sim

    def to_file(self, output_filename: str) -> None:
        """Write an object to disk at the filename given.

        Parameters:
            output_filename (str):The filename to write to.  Files get overwritten
                                  by default.

        Returns:
            None
        """

        with open(output_filename, 'wb') as fout:
            pickle.dump(self, fout, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def from_file(input_filename: str) -> SimulationFactory:
        """Read an object from disk from the filename given.

        Parameters:
            input_filename (str):The filename to read from.

        Returns:
            obj(Object):The python object read from disk
        """

        log.info('Reading data from %s...', input_filename)
        with open(input_filename, 'rb') as fin:
            payload = pickle.load(fin)

        return payload
