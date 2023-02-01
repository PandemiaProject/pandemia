from pandemia.world.region import Region, VectorRegion
from collections import namedtuple
import pytest

def test_vector_region_constructor():
    Args = namedtuple("Args", [
            "id",
            "name",
            "ticks_in_week",
            "number_of_activities",
            "number_of_agents",
            "number_of_locations",
            "max_num_activity_locations"
    ])

    passing_cases = [
        Args(
            id=1,
            name="a",
            ticks_in_week=21,
            number_of_activities=3,
            number_of_agents=4,
            number_of_locations=4,
            max_num_activity_locations=12
        ),
        Args(
            id=1,
            name="a",
            ticks_in_week=21,
            number_of_activities=3,
            number_of_agents=4,
            number_of_locations=4,
            max_num_activity_locations=12
        )
    ]

    failing_cases = [
        Args(
            id=1,
            name="a",
            ticks_in_week=-5, # Negative number of ticks
            number_of_activities=3,
            number_of_agents=4,
            number_of_locations=4,
            max_num_activity_locations=12
        ),
        Args(
            id=1,
            name="a",
            ticks_in_week=21,
            number_of_activities=3,
            number_of_agents=4,
            number_of_locations=4,
            # Implausible large value given
            # `number_of_activities=3` and
            # `number_of_agents=4``
            max_num_activity_locations=120000
        ),
    ]

    for arg in passing_cases:
        VectorRegion(
            id=arg.id,
            name=arg.name,
            ticks_in_week=arg.ticks_in_week,
            number_of_activities=arg.number_of_activities,
            number_of_agents=arg.number_of_agents,
            number_of_locations=arg.number_of_locations,
            max_num_activity_locations=arg.max_num_activity_locations
        )
        assert True

    for arg in failing_cases:
        with pytest.raises(ValueError, match="(Out of range value|negative dimensions are not allowed)"):
            VectorRegion(
                id=arg.id,
                name=arg.name,
                ticks_in_week=arg.ticks_in_week,
                number_of_activities=arg.number_of_activities,
                number_of_agents=arg.number_of_agents,
                number_of_locations=arg.number_of_locations,
                max_num_activity_locations=arg.max_num_activity_locations
            )
