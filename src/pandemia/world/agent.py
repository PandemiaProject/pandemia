
import logging
import uuid

from pandemia.world.location import Location

log = logging.getLogger("agent")

class Agent:
    """Represents a location, for example an area of land or a building.

    Parameters
    ----------
    age : `int`
        The age of the individual, in years.
    """

    uuid = None
    """A universally unique identifier for this agent (`str`).
    """

    age = None
    """The age of the individual, in years (`int`).
    """

    weekly_routine = None
    """A sequence of activites, representing a weekly routine (`list[str]`). This specifies which
    activities this agent performs and in which order. For example, suppose this agent performs two
    activities, "Home" and "Work". Then their weekly routine might be:

            ["Home", "Home", "Home", "Home", "Work", "Home", "Home", "Work", "Home",
            "Home", "Work", "Home", "Home", "Work", "Home", "Home", "Work", "Home",
            "Home", "Home", "Home"]

    In this example, each string refers to an eight-hour time interval, since it will be assumed
    that the routine divides the week into time intervals of equal length.
    """

    activity_locations = None
    """For each string appearing in the weekly routine, a list of objects of type Location,
    indicating at which locations it is possible for the agent to perform that activity
    (`dict[str, list[Location]]`). In the weekly_routine example, the agent would require two lists,
    one for activity "Home" and one for activity "Work". The list for activity "Home" might consist
    of a single Location of an appropriate type, for example "House".
    """

    activity_location_weights = None
    """For each of the activity_locations lists, an associated list of weights, indicating how
    likely it is that the agent will choose the corresponding locations when performing each
    activity (`dict[str, list[float]]`). These lengths of these lists of weights should therefore be
    the same as the lengths of the lists of locations, for each activity. The weights need not sum
    to one."""


    def __init__(self, age):

        self.uuid: str = uuid.uuid4().hex
        self.age: int = age
        self.weekly_routine: list[str] = []
        self.activity_locations: dict[str, list[Location]] = {}
        self.activity_location_weights: dict[str, list[float]] = {}
