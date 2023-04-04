"""Builds a Simulator object. The SimulationFactory object builds a Simulator object gradually, by
adding the various model layers and components."""

# Allows classes to return their own type, e.g. from_file below
from __future__ import annotations

import logging
import uuid
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
from .components.policy_maker_model import PolicyMakerModel

log = logging.getLogger("sim_state")

class SimulationFactory:
    """Class that allows for gradual composition of a number of components, eventually outputting
    a Simulator object that can be used to run simulations with the given config."""

    def __init__(self, config: Config, vector_world: VectorWorld):

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
        self.vector_world = vector_world

        # Components
        self.policy_maker_model = None
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

    def set_policy_maker_model(self, policy_maker_model: PolicyMakerModel) -> None:
        """Sets seasonal effects model"""
        self.policy_maker_model = policy_maker_model

    def build_clock(self):
        """Builds the simulation clock"""

        config = self.config

        # Create clock
        _clock = Clock(config['tick_length_s'], config['simulation_length_days'], config['epoch'])
        log.info("New clock created at %s, tick_length = %i, simulation_days = %i",
                 _clock.epoch, _clock.tick_length_s, _clock.simulation_length_days)
        self.set_clock(_clock)

        # Check that time discretization in world (weekly routines) matches that of clock
        ticks_in_week_clock = self.clock.ticks_in_week
        for vector_region in self.vector_world.vector_regions:
            ticks_in_week_region = vector_region.weekly_routines.shape[1]
            assert ticks_in_week_clock == ticks_in_week_region

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
        health_model = instantiate_class("pandemia.components.health_model",
                                        health_model_class, health_model_config, scale_factor,
                                        self.clock)
        self.set_health_model(health_model)

        # Create movement model
        movement_model_class = config['movement_model.__type__']
        movement_model_config = config.subconfig('movement_model')
        movement_model = instantiate_class("pandemia.components.movement_model",
                                        movement_model_class, movement_model_config)
        self.set_movement_model(movement_model)

        # Create hospitalization and death model
        hospitalization_and_death_model_class = config['hospitalization_and_death_model.__type__']
        hospitalization_and_death_model_config = config.subconfig('hospitalization_and_death_model')
        hospitalization_and_death_model =\
            instantiate_class("pandemia.components.hospitalization_and_death_model",
                            hospitalization_and_death_model_class,
                            hospitalization_and_death_model_config)
        self.set_hospitalization_and_death_model(hospitalization_and_death_model)

        # Create testing and contact tracing model
        testing_and_contact_tracing_model_class =\
            config['testing_and_contact_tracing_model.__type__']
        testing_and_contact_tracing_model_config =\
            config.subconfig('testing_and_contact_tracing_model')
        testing_and_contact_tracing_model =\
            instantiate_class("pandemia.components.testing_and_contact_tracing_model",
                            testing_and_contact_tracing_model_class,
                            testing_and_contact_tracing_model_config)
        self.set_testing_and_contact_tracing_model(testing_and_contact_tracing_model)

        # Create vaccination model
        vaccination_model_class = config['vaccination_model.__type__']
        vaccination_model_config = config.subconfig('vaccination_model')
        vaccination_model = instantiate_class("pandemia.components.vaccination_model",
                                            vaccination_model_class, vaccination_model_config,
                                            self.clock,
                                            health_model.number_of_strains,
                                            health_model.number_of_rho_immunity_outcomes,
                                            health_model.immunity_length,
                                            health_model.immunity_period_ticks)
        self.set_vaccination_model(vaccination_model)

        # Create regional mixing model
        travel_model_class = config['travel_model.__type__']
        travel_config = config.subconfig('travel_model')
        travel = instantiate_class("pandemia.components.travel_model",
                                            travel_model_class, travel_config,
                                            scale_factor,
                                            health_model.number_of_strains,
                                            self.vector_world)
        self.set_travel_model(travel)

        # Create policy maker model
        policy_maker_class = config['policy_maker_model.__type__']
        policy_maker_config = config.subconfig('policy_maker_model')
        policy_maker = instantiate_class("pandemia.components.policy_maker_model",
                                policy_maker_class, policy_maker_config, scale_factor,
                                self.clock,
                                self.vector_world.number_of_regions,
                                vaccination_model.number_of_vaccines,
                                vaccination_model.age_groups)
        self.set_policy_maker_model(policy_maker)

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
                        self.policy_maker_model,
                        telemetry_bus)

        return sim
