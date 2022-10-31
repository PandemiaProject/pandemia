"""World factories build a World object.

Worlds consist of locations, agents and activities.

They provide a set of data for the simulation to run upon.
"""

from pandemia.world import World

class WorldFactory:
    """Generic world factory that outputs a World object."""

    def get_world(self) -> World:
        """Gets world"""
        pass
