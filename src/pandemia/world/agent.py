"""Represents an individual, assigning to each an age, weekly routine and for each activity a list
of locations at which they can perform that activity"""

import logging
import uuid

from pandemia.world.location import Location

log = logging.getLogger("agent")

class Agent:
    """Represents an individual"""

    def __init__(self, age):

        self.uuid = uuid.uuid4().hex
        self.age: int = age
        self.weekly_routine: list[str] = []
        self.activity_locations: dict[str, list[Location]] = {}
        self.activity_location_weights: dict[str, list[float]] = {}
