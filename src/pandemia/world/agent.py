"""Represents a individual"""

import logging
import uuid

from .location import Location

log = logging.getLogger("agent")

class Agent:

    def __init__(self, age):
        """Represents a location, for example an area of land or a building.

        Attributes:
           uuid (str): A universally unique identifier for this agent.
           age (int): The age of the individual, in years.
           weekly_routine (list[str]): A list of activites, representing a weekly routine. This
              specifies which activities this agent performs and in which order. For example,
              suppose this agent performs two activities, "Home" and "Work". Then their weekly
              routine might be:

                  ["Home", "Home", "Home", "Home", "Work", "Home", "Home", "Work", "Home",
                   "Home", "Work", "Home", "Home", "Work", "Home", "Home", "Work", "Home",
                   "Home", "Home", "Home"]

              In this example, it is to be understood that each string refers to an eight-hour time
              interval, since it is to be assumed that the week is divided into time intervals of
              equal length.
          activity_locations (dict[str, list[Location]]): For each string appearing in their weekly
              routine, a list of objects of type Location is assigned, indicating at which locations
              it is posisble for the agent to perform that activity. In the above example, an agent
              would require two lists, one for activity "Home" and one for activity "Work". The list
              for activity "Home" might consist of a single Location of an appropriate type, for
              example "House".
          activity_location_weights (dict[str, list[float]]): For each of the above lists, an
              associated list of weights, indicating how likely it is that the agent will choose the
              corresponding locations when performing each activity. These lists of weights should
              therefore be the same length as the lists of locations, for each activity. The weights
              need not sum to 1.
        """

        self.uuid: str = uuid.uuid4().hex
        self.age: int = age
        self.weekly_routine: list[str] = []
        self.activity_locations: dict[str, list[Location]] = {}
        self.activity_location_weights: dict[str, list[float]] = {}
