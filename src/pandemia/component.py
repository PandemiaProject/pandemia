"""Components for the simulation"""

import logging
from typing import Optional
from abc import ABC, abstractmethod

from .config import Config
from .messagebus import MessageBus
from ctypes import cdll
import platform
import os
ext=".dll" if platform.system() == 'Windows' else ".so"
log = logging.getLogger("component")

#pylint: disable=attribute-defined-outside-init
class Component(ABC):
    """A simulation component. Submodels, for example of movement and health, are represented as
    objects of this class."""

    def __init__(self, component_config: Config):
        """Initialize the component"""

        self.lib = cdll.LoadLibrary(os.path.split(__file__)[0]+"/C/build/default_model_functions"+ext)
        self.config = component_config
        self.telemetry_bus: Optional[MessageBus] = None

    def set_telemetry_bus(self, telemetry_bus: Optional[MessageBus]) -> None:
        """Sets telemetry bus for this component"""
        self.telemetry_bus = telemetry_bus

    def report(self, topic, *args, **kwargs):
        """Publish a message to the telemetry bus, informing reporters of some interesting event
        or statistic."""

        if self.telemetry_bus is None:
            return

        self.telemetry_bus.publish(topic, *args, **kwargs)

    @abstractmethod
    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""
        pass

    @abstractmethod
    def initial_conditions(self, vector_region):
        """Define the initial state of the vector_regions"""
        pass

