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
    vaccination_model
)

"""
from pandemia.components.health_model.default_health_model import DefaultHealthModel
from pandemia.components.health_model import HealthModel
from pandemia.components.hospitalization_and_death_model import HospitalizationAndDeathModel
from pandemia.components.health_model.void_health_model import VoidHealthModel
from pandemia.components.hospitalization_and_death_model.default_hospitalization_and_death_model import DefaultHospitalizationAndDeathModel
from pandemia.components.hospitalization_and_death_model.void_hospitalization_and_death_model import VoidHospitalizationAndDeathModel
from pandemia.components.movement_model import MovementModel
from pandemia.components.movement_model.default_movement_model import DefaultMovementModel
from pandemia.components.movement_model.void_movement_model import VoidMovementModel
from pandemia.components.policy_maker_model import PolicyMakerModel
from pandemia.components.policy_maker_model.default_policy_maker_model import DefaultPolicyMakerModel
from pandemia.components.policy_maker_model.optimization_policy_maker_model import OptimizationPolicyMakerModel
from pandemia.components.policy_maker_model.random_policy_maker_model import RandomPolicyMakerModel
from pandemia.components.policy_maker_model.validation_policy_maker_model import ValidationPolicyMakerModel
from pandemia.components.policy_maker_model.void_policy_maker_model import VoidPolicyMakerModel
from pandemia.components.seasonal_effects_model import SeasonalEffectsModel
from pandemia.components.seasonal_effects_model.default_seasonal_effects_model import DefaultSeasonalEffectsModel
from pandemia.components.seasonal_effects_model.void_seasonal_effects_model import VoidSeasonalEffectsModel
from pandemia.components.testing_and_contact_tracing_model import TestingAndContactTracingModel
from pandemia.components.testing_and_contact_tracing_model.default_testing_and_contact_tracing_model import DefaultTestingAndContactTracingModel
from pandemia.components.testing_and_contact_tracing_model.void_testing_and_contact_tracing_model import VoidTestingAndContactTracingModel
from pandemia.components.travel_model import TravelModel
from pandemia.components.travel_model.default_travel_model import DefaultTravelModel
from pandemia.components.travel_model.void_travel_model import VoidTravelModel
from pandemia.components.vaccination_model import VaccinationModel
from pandemia.components.vaccination_model.default_vaccination_model import DefaultVaccinationModel
from pandemia.components.vaccination_model.void_vaccination_model import VoidVaccinationModel
"""



from pandemia.config import Config
from pandemia.world import region
import pytest
from pkgutil import walk_packages
import importlib
from icecream import ic

@pytest.fixture
def list_of_abstract_components():
    """These classes _should_ be abstract"""
    return [
        Component,
        health_model.HealthModel,
        hospitalization_and_death_model.HospitalizationAndDeathModel,
        movement_model.MovementModel,
        policy_maker_model.PolicyMakerModel,
        seasonal_effects_model.SeasonalEffectsModel,
        testing_and_contact_tracing_model.TestingAndContactTracingModel,
        travel_model.TravelModel,
        vaccination_model.VaccinationModel
    ]

def _all_subclasses(cls):
    """
    Note that this will only find subclasses that have already been imported
    See https://stackoverflow.com/a/3862957/3837936
    """
    return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in _all_subclasses(c)])

@pytest.fixture
def all_components():
    """
    This will return `Component`, and all known subclasses (both abstract and concrete)
    """
    # Find and import all the subpackages
    search_path = [pandemia.__path__[0] + "/components"]

    for _, module_name, _ in walk_packages(search_path, pandemia.__name__ + '.components.'):
        print(module_name)
        importlib.import_module(module_name)

    # Now search for subclasses of `Component`
    subclasses = _all_subclasses(Component)
    return subclasses



@pytest.fixture
def sample_config():

    return Config(_dict={})
@pytest.fixture

def sample_vector_region():
    def v_region():
        return region.VectorRegion(
            id=1,
            name="test_vector_region",
            ticks_in_week=21,
            number_of_activities=2,
            number_of_agents=2,
            number_of_locations=1,
            max_num_activity_locations=2
        )

    return v_region

def test_abstract_components(list_of_abstract_components, sample_config):
    """Test that the top level components are abstract and cannot be instantiated"""
    for comp in list_of_abstract_components:
        print(str(comp))
        with pytest.raises(TypeError, match="Can't instantiate abstract class"):
            comp(sample_config)

def test_component_constructor(all_components, sample_config):
    for comp in all_components:
        print(comp)
        comp(sample_config)

    assert False

def test_vectorize_component(all_components, sample_config, sample_vector_region):
    for Comp in all_components:
        comp = Comp(sample_config)
        v_region = sample_vector_region()
        before = len(dir(v_region))
        comp.vectorize_component(v_region)
        after = len(dir(v_region))
        ic(after - before, Comp)

    assert False

@pytest.mark.skip("Test not implemented yet")
def test_initial_conditions():
    assert False

@pytest.mark.skip("Test not implemented yet")
def test_component_dynamics():
    assert False
