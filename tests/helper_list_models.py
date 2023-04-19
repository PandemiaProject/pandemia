from pandemia.components.health_model.default_health_model import DefaultHealthModel
from pandemia.components.health_model.safir_health_model import SafirHealthModel
from pandemia.components.health_model.void_health_model import VoidHealthModel
from pandemia.components.hospitalization_and_death_model.default_hospitalization_and_death_model import (
    DefaultHospitalizationAndDeathModel,
)
from pandemia.components.hospitalization_and_death_model.void_hospitalization_and_death_model import (
    VoidHospitalizationAndDeathModel,
)
from pandemia.components.movement_model.default_movement_model import (
    DefaultMovementModel,
)
from pandemia.components.movement_model.void_movement_model import VoidMovementModel
from pandemia.components.policy_maker_model.default_policy_maker_model import (
    DefaultPolicyMakerModel,
)
from pandemia.components.policy_maker_model.optimization_policy_maker_model import (
    OptimizationPolicyMakerModel,
)
from pandemia.components.policy_maker_model.random_policy_maker_model import (
    RandomPolicyMakerModel,
)
from pandemia.components.policy_maker_model.validation_policy_maker_model import (
    ValidationPolicyMakerModel,
)
from pandemia.components.policy_maker_model.void_policy_maker_model import (
    VoidPolicyMakerModel,
)
from pandemia.components.seasonal_effects_model.default_seasonal_effects_model import (
    DefaultSeasonalEffectsModel,
)
from pandemia.components.seasonal_effects_model.void_seasonal_effects_model import (
    VoidSeasonalEffectsModel,
)
from pandemia.components.testing_and_contact_tracing_model.default_testing_and_contact_tracing_model import (
    DefaultTestingAndContactTracingModel,
)
from pandemia.components.testing_and_contact_tracing_model.void_testing_and_contact_tracing_model import (
    VoidTestingAndContactTracingModel,
)
from pandemia.components.travel_model.default_travel_model import DefaultTravelModel
from pandemia.components.travel_model.safir_travel_model import SafirTravelModel
from pandemia.components.travel_model.void_travel_model import VoidTravelModel
from pandemia.components.vaccination_model.default_vaccination_model import (
    DefaultVaccinationModel,
)
from pandemia.components.vaccination_model.safir_vaccination_model import (
    SafirVaccinationModel,
)
from pandemia.components.vaccination_model.void_vaccination_model import (
    VoidVaccinationModel,
)

from pandemia.config import Config
from pandemia.clock import Clock
from pandemia.world.world_factory.control_world_factory import ControlWorldFactory

default_conf = Config("Scenarios/Test/test_all_components.yaml")
void_conf = Config(_dict={})

validation_conf = Config("Scenarios/Validation/heterogeneous_validation_four_strains.yaml")
safir_conf = Config("Scenarios/Safir/safir.yaml")

world_factory_conf = Config(
    _dict={
        "random_seed": 1,
        # Time discretization
        "ticks_in_day_world": 144,
        # The order in which the region names are listed here determines the index of that region referred
        # to in vectors and matrices appearing in other components:
        "region_names": [
            "Test Region 0",
            "Test Region 1",
            "Test Region 2",
            "Test Region 3",
        ],
        # Locations
        "location_type_counts": {
            "House": 25,
            "Hospital": 2,
            "Cemetery": 1,
            "Other_Location": 22,
        },
        # Agents
        "number_of_agents_per_region": 10,
        "max_age": 100,
        # Activities
        "activities": ["Home_Activity", "Other_Activity"],
        "activity_location_types": {
            "Home_Activity": ["House"],
            "Other_Activity": ["Other_Location"],
        },
        "number_of_locations_for_agent_by_activity": {
            "Home_Activity": 1,
            "Other_Activity": 10,
        },
        "contacts": 5,
        "contact_hours": 8,
        "travel_prob_per_day": 0.0001,
    }
)


def get_clock():
    return Clock(60 * 60 * 8, 2, "1st January 2020")


def get_vector_world():
    world_factory = ControlWorldFactory(
        config=world_factory_conf,
        scale_factor=1
    )

    return world_factory.get_world().vectorize_world()


all_model_list = [
    DefaultHealthModel(
        config=default_conf.subconfig("health_model"), scale_factor=1, clock=get_clock()
    ),
    SafirHealthModel(
        config=safir_conf.subconfig("health_model"), scale_factor=1, clock=get_clock()
    ),
    VoidHealthModel(config=void_conf, scale_factor=1, clock=get_clock()),
    DefaultHospitalizationAndDeathModel(
        config=default_conf.subconfig("hospitalization_and_death_model")
    ),
    VoidHospitalizationAndDeathModel(config=void_conf),
    DefaultMovementModel(config=default_conf.subconfig("movement_model")),
    VoidMovementModel(config=void_conf),
    DefaultPolicyMakerModel(
        config=Config(_dict={"policy_data_filepath": "some_file"}),
        scale_factor=1,
        clock=get_clock(),
        number_of_regions=4,
        number_of_vaccines=2,
        age_groups=[0, 18, 65],
    ),
    OptimizationPolicyMakerModel(
        config=void_conf,
        scale_factor=1,
        clock=get_clock(),
        number_of_regions=4,
        number_of_vaccines=2,
        age_groups=[0, 18, 65],
    ),
    RandomPolicyMakerModel(
        config=default_conf.subconfig("policy_maker_model"),
        scale_factor=1,
        clock=get_clock(),
        number_of_regions=4,
        number_of_vaccines=2,
        age_groups=[0, 18, 65],
    ),
    ValidationPolicyMakerModel(
        config=validation_conf.subconfig("policy_maker_model"),
        scale_factor=1,
        clock=get_clock(),
        number_of_regions=4,
        number_of_vaccines=2,
        age_groups=[0, 18, 65],
    ),
    VoidPolicyMakerModel(
        config=void_conf,
        scale_factor=1,
        clock=get_clock(),
        number_of_regions=4,
        number_of_vaccines=2,
        age_groups=[0, 18, 65],
    ),
    DefaultSeasonalEffectsModel(
        config=default_conf.subconfig("seasonal_effects_model"),
        vector_world=get_vector_world(),
        clock=get_clock(),
    ),
    VoidSeasonalEffectsModel(
        config=void_conf,
        vector_world=get_vector_world(),
        clock=get_clock(),
    ),
    DefaultTestingAndContactTracingModel(
        config=default_conf.subconfig("testing_and_contact_tracing_model")
    ),
    VoidTestingAndContactTracingModel(config=void_conf),
    DefaultTravelModel(
        config=default_conf.subconfig("travel_model"),
        scale_factor=1,
        number_of_strains=2,
        vector_world=get_vector_world(),
    ),
    SafirTravelModel(
        config=safir_conf.subconfig("travel_model"),
        scale_factor=1,
        number_of_strains=2,
        vector_world=get_vector_world(),
    ),
    VoidTravelModel(
        config=void_conf, scale_factor=1, number_of_strains=2, number_of_regions=4
    ),
    DefaultVaccinationModel(
        config=default_conf.subconfig("vaccination_model"),
        clock=get_clock(),
        number_of_strains=2,
        number_of_rho_immunity_outcomes=2,
        immunity_length=14,
        immunity_period_ticks=14 * 3,
    ),
    SafirVaccinationModel(
        config=safir_conf.subconfig("vaccination_model"),
        clock=get_clock(),
        number_of_strains=2,
        number_of_rho_immunity_outcomes=2,
        immunity_length=14,
        immunity_period_ticks=14 * 3,
    ),
    VoidVaccinationModel(
        config=void_conf,
        clock=get_clock(),
        number_of_strains=2,
        number_of_rho_immunity_outcomes=2,
        immunity_length=14,
        immunity_period_ticks=14 * 3,
    ),
]
