"""Components for the simulation"""

import logging
from typing import Optional

from pandemia.config import Config
from pandemia.messagebus import MessageBus

log = logging.getLogger("component")

#pylint: disable=attribute-defined-outside-init
class Component:
    """A simulation component. Submodels, for example of movement and health, are represented as
    objects of this class."""

    def __init__(self, component_config: Config):
        """Initialize the component"""

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
