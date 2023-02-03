import pandemia
from pandemia.component import Component
from pandemia.components import (
    health_model,
    hospitalization_and_death_model,
    movement_model,
    policy_maker_model,
    seasonal_effects_model,
    testing_and_contact_tracing_model,
    travel_model,
    vaccination_model,
)

from pandemia.random_tools import Random


from pandemia.config import Config
from pandemia.world import region
import pytest
from pkgutil import walk_packages
import importlib
from icecream import ic

from helper_list_models import all_model_list

abstract_classes = [
    Component,
    health_model.HealthModel,
    hospitalization_and_death_model.HospitalizationAndDeathModel,
    movement_model.MovementModel,
    policy_maker_model.PolicyMakerModel,
    seasonal_effects_model.SeasonalEffectsModel,
    testing_and_contact_tracing_model.TestingAndContactTracingModel,
    travel_model.TravelModel,
    vaccination_model.VaccinationModel,
]


@pytest.fixture(
    scope="module",
    params=all_model_list,
    ids=[type(obj).__name__ for obj in all_model_list],
)
def concrete_model(request):
    return request.param


def _recurse_subclasses(cls):
    """
    Note that this will only find subclasses that have already been imported
    See https://stackoverflow.com/a/3862957/3837936
    """
    return set(cls.__subclasses__()).union(
        [s for c in cls.__subclasses__() for s in _recurse_subclasses(c)]
    )


@pytest.fixture
def find_subclasses():
    """
    This will return `Component`, and all known subclasses (both abstract and concrete)
    """
    # Find and import all the subpackages
    search_path = [pandemia.__path__[0] + "/components"]

    for _, module_name, _ in walk_packages(
        search_path, pandemia.__name__ + ".components."
    ):
        print(module_name)
        importlib.import_module(module_name)

    # Now search for subclasses of `Component`
    subclasses = _recurse_subclasses(Component)
    return subclasses


@pytest.fixture
def sample_vector_region():
    """Create a new instance of VectorRegion for each call"""

    def v_region():
        vr = region.VectorRegion(
            id=1,
            name="test_vector_region",
            ticks_in_week=21,
            number_of_activities=2,
            number_of_agents=2,
            number_of_locations=1,
            max_num_activity_locations=2,
        )
        vr.prng = Random(123)
        return vr

    return v_region


@pytest.mark.parametrize(
    "class_name",
    abstract_classes,
)
def test_abstract_components(class_name):
    """Test that the top level components are abstract and cannot be instantiated. These classes _should_ be abstract"""
    with pytest.raises(TypeError, match="Can't instantiate abstract class"):
        class_name(Config(_dict={}))


def test_sample_exist_for_each_concrete_class(find_subclasses):
    """
    This test asserts that the list `helper_list_models.all_model_list` contains exactly one of each
    available concrete subclass of `Component`. If a new subclass is created and an example not added to
    `helper_list_models.all_model_list` then this test will fail.
    """
    actual_subclasses = set(find_subclasses) - set(abstract_classes)
    expected_subclasses = set([type(obj) for obj in all_model_list])

    diff = actual_subclasses.symmetric_difference(expected_subclasses)
    assert len(diff) == 0


def test_vectorize_component(concrete_model, sample_vector_region):
    """
    For all Component subclasses - test calling the `vectorize_component` method.
    In a simulation this is typically the first method called, after construction.

    **At present the only positive assertion that this has worked as expected is the change in the number of attributes on the `vector_region`.**

    NOTE: The expected values have _only_ been obtained by running this test. The have not been checked by manually inspecting the code.
    """

    # Find expected failure cases
    xfail_if = {
        "DefaultMovementModel": 87,
        "ValidationPolicyMakerModel": 89,
    }

    expected_fail_issue = xfail_if.get(type(concrete_model).__name__, None)

    if expected_fail_issue:
        expected_fail_reason = f"Expected fail. See issue https://github.com/PandemiaProject/pandemia/issues/{expected_fail_issue}"
        pytest.xfail(reason=expected_fail_reason)

    # Now do the actual tests for the remaining cases
    all_excepted_values = {
        "DefaultHealthModel": 27,
        "VoidHealthModel": 9,
        "DefaultHospitalizationAndDeathModel": 6,
        "VoidHospitalizationAndDeathModel": 0,
        "DefaultMovementModel": 1,
        "VoidMovementModel": 9,
        "DefaultPolicyMakerModel": 0,
        "OptimizationPolicyMakerModel": 0,
        "RandomPolicyMakerModel": 0,
        "ValidationPolicyMakerModel": 1,
        "VoidPolicyMakerModel": 0,
        "DefaultSeasonalEffectsModel": 2,
        "VoidSeasonalEffectsModel": 1,
        "DefaultTestingAndContactTracingModel": 8,
        "VoidTestingAndContactTracingModel": 0,
        "DefaultTravelModel": 2,
        "VoidTravelModel": 1,
        "DefaultVaccinationModel": 5,
        "VoidVaccinationModel": 0,
    }

    expected_extra_attributes = all_excepted_values[type(concrete_model).__name__]

    v_region = sample_vector_region()
    before = len(dir(v_region))
    concrete_model.vectorize_component(v_region)
    after = len(dir(v_region))
    actual_extra_attributes = after - before
    ic(actual_extra_attributes, concrete_model)
    assert actual_extra_attributes == expected_extra_attributes


def test_initial_conditions(concrete_model, sample_vector_region):
    """
    For all Component subclasses - test calling the `initial_conditions` method.

    The `vectorize_component` method is called first, has it is typically the first method called, after construction.
    Any errors from `vectorize_component` are ignored in this test.

    **At present the NO positive assertion that this has worked as expected. Only the absence of an error in
    `initial_conditions` is sufficient for this test to pass.
    """

    # Find expected failure cases
    xfail_if = {
        "DefaultHealthModel": 86,
        "DefaultMovementModel": 88,
        "VoidMovementModel": 88,
        "DefaultPolicyMakerModel": 89,
        "OptimizationPolicyMakerModel": 89,
        "RandomPolicyMakerModel": 89,
        "ValidationPolicyMakerModel": 89,
        "DefaultTestingAndContactTracingModel": 90,
        "DefaultTravelModel": 90,
    }

    expected_fail_issue = xfail_if.get(type(concrete_model).__name__, None)

    if expected_fail_issue:
        expected_fail_reason = f"Expected fail. See issue https://github.com/PandemiaProject/pandemia/issues/{expected_fail_issue}"
        pytest.xfail(reason=expected_fail_reason)

    # Now do the actual tests for the remaining cases
    v_region = sample_vector_region()

    # Ignore errors in `vectorize_component`
    try:
        concrete_model.vectorize_component(v_region)
    except Exception:
        pass

    concrete_model.initial_conditions(v_region)
