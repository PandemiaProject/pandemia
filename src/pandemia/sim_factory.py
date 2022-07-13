"""Represents the simulation state.  Is built gradually by the various model stages, and then
ingested by the simulator as it runs."""

# Allows classes to return their own type, e.g. from_file below
from __future__ import annotations

import logging
import uuid
import pickle
from datetime import datetime

from pandemia.config import Config
from pandemia.clock import Clock
from pandemia.version import VERSION
from pandemia.world import VectorWorld
from pandemia.simulator import Simulator
from pandemia.components.seasonal_effects_model import SeasonalEffectsModel
from pandemia.components.health_model import HealthModel
from pandemia.components.movement_model import MovementModel
from pandemia.components.hospitalization_and_death_model import HospitalizationAndDeathModel
from pandemia.components.testing_and_contact_tracing_model import TestingAndContactTracingModel
from pandemia.components.vaccination_model import VaccinationModel
from pandemia.components.regional_mixing_model import RegionalMixingModel
from pandemia.components.input_model import InputModel

log = logging.getLogger("sim_state")

class SimulationFactory:
    """Class that allows for gradual composition of a number of components, eventually outputting
    a simulator object that can be used to run simulations with the config given."""

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
        self.regional_mixing_model = None

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

    def set_regional_mixing_model(self, regional_mixing_model: RegionalMixingModel) -> None:
        """Sets regional mixing model"""
        self.regional_mixing_model = regional_mixing_model

    def set_input_model(self, input_model: InputModel) -> None:
        """Sets seasonal effects model"""
        self.input_model = input_model

    def new_sim(self, telemetry_bus):
        """Return a new simulator based on the config above.

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
                        self.regional_mixing_model,
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
