import sys
import pandas as pd
import pandemia
from pandemia.__main__ import main as pandemia_main
import pytest
import yaml
from pathlib import Path
from typing import Union

from icecream import ic


def _run_end_to_end_simulation(
    input_scenario_config, expected_results_csv, output_dir: Path, output_fname: str
):
    """
    This function will run a complete simulation. Most of the end to end test functions are a wrapper for
    this function.

    Because these tests are slow to run, it is possible that developer might wish to run tests separately.
    (If these tests run quickly (eg <10sec) then it might be more appropriate to loop through all of the
    possible scenarios in the single test function).

    :params:
    @input_scenario_config: The path to the scenario file. This can be absolute path, or relative to the
    root of the repo.
    @expected_results_csv: The path to a csv file containing the expected results. This can be absolute
    path, or relative to the the root of the repo. The test passes or fails based on the comparison of
    this csv file and the output of the simulation will be with this csv, using the
    `pandas.DataFrame.compare` function.
    @output_dir A `pathlib.Path` object, representing a directory to write outputs.
    Typically this the `tmp_path` fixture (https://docs.pytest.org/en/latest/how-to/tmp_path.html#the-tmp-path-fixture)
    @output_fname the output filename. If the value does not end with the suffix ".csv" it will be appended automatically.

    """
    fname_altered_scenario = str(output_dir / "test_altered_scenario.yaml")

    if output_fname.endswith(".csv"):
        fname_actual_results = str(output_dir / output_fname)
    else:
        fname_actual_results = str(output_dir / f"{output_fname}.csv")

    ic(fname_altered_scenario)
    ic(fname_actual_results)

    output_location = {
        "reporters": {"csv.StrainCounts": {"filename": fname_actual_results}}
    }

    scenario = _update_config(input_scenario_config, output_location)

    with open(fname_altered_scenario, "w") as out_yaml:
        yaml.dump(scenario, out_yaml)

    # Read the expected results
    # Do this before the main simulation, so that the test fails quickly in the case of a
    # missing/unreadable `expected_results_csv` file.
    expected_df = pd.read_csv(expected_results_csv)

    # Fake the commandline args and call Pandemia.main()
    sys.argv = [pandemia.__file__, fname_altered_scenario]
    # pandemia.__main__.main()
    pandemia_main()

    # Actual results
    actual_df = pd.read_csv(fname_actual_results)

    # Check for differences in the results
    diff_df = expected_df.compare(actual_df)

    assert len(diff_df) == 0
    assert len(diff_df.index) == 0
    assert len(diff_df.columns) == 0


def _update_config(input_scenario: Union[dict, str, Path], delta_config: dict) -> dict:
    """
    @input_scenario either a dict or a path representing a input configuration. If it
    @delta_config A dict which a set of changes to apply to input_scenario
    @returns A config dict, which has been updated using delta_confg
    """

    if isinstance(input_scenario, dict):
        scenario = input_scenario
    else:
        # Update the config with temporary dir for the output files
        with open(input_scenario, "r") as in_yaml:
            scenario = yaml.safe_load(in_yaml)

    scenario.update(delta_config)

    return scenario


@pytest.mark.slow
def test_end_to_end_global(tmp_path, request):
    input_scenario = "Scenarios/Test/test_global_config.yaml"
    # expect_results = "tests/e2e_expected_outputs/strain_counts.csv"
    expect_results = "tests/e2e_expected_outputs/strain_counts_test_global_config.csv"

    _run_end_to_end_simulation(
        input_scenario, expect_results, tmp_path, request.node.originalname
    )


@pytest.mark.known_failing
def test_end_to_end_all_components(tmp_path, request):
    input_scenario = "Scenarios/Test/test_all_components.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts_test_all_components.csv"

    _run_end_to_end_simulation(
        input_scenario, expect_results, tmp_path, request.node.originalname
    )


@pytest.mark.slow
def test_end_to_end_all_void(tmp_path, request):
    input_scenario = "Scenarios/Test/test_void_all.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts_test_void_all.csv"

    _run_end_to_end_simulation(
        input_scenario, expect_results, tmp_path, request.node.originalname
    )


@pytest.mark.slow
def test_e2e_health_model(tmp_path, request):
    input_scenario = "Scenarios/Test/test_e2e_health_model.yaml"
    expect_results = (
        "tests/e2e_expected_outputs/test_e2e_health_model.csv"
    )

    _run_end_to_end_simulation(
        input_scenario, expect_results, tmp_path, request.node.originalname
    )


# @pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_hospitalization_and_death_model(tmp_path, request):
    # input_scenario = "Scenarios/Test/test_e2e_hospitalization_and_death_model.yaml"
    input_scenario = "Scenarios/Test/test_e2e_health_model.yaml"
    expect_results = "tests/e2e_expected_outputs/test_e2e_hospitalization_and_death_model.csv"

    delta_config = {
        "hospitalization_and_death_model": {
            "__type__": "default_hospitalization_and_death_model.DefaultHospitalizationAndDeathModel",
            "hospital_threshold": 0.5,
            "hospital_location_type": "Hospital",
            "cemetery_location_type": "Cemetery",
        }
    }

    scenario = _update_config(input_scenario, delta_config)
    _run_end_to_end_simulation(
        scenario, expect_results, tmp_path, request.node.originalname
    )


# @pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_input_model(tmp_path, request):
    # input_scenario = "Scenarios/Test/test_e2e_input_model.yaml"
    input_scenario = "Scenarios/Test/test_e2e_health_model.yaml"
    expect_results = "tests/e2e_expected_outputs/test_e2e_input_model.csv"

    delta_config = {
        "input_model": {"__type__": "default_input_model.DefaultInputModel"}
    }

    # run_end_to_end_simulation(input_scenario, expect_results, tmp_path, request.node.originalname)
    scenario = _update_config(input_scenario, delta_config)
    _run_end_to_end_simulation(
        scenario, expect_results, tmp_path, request.node.originalname
    )


# @pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_movement_model(tmp_path, request):
    # input_scenario = "Scenarios/Test/test_e2e_movement_model.yaml"
    input_scenario = "Scenarios/Test/test_e2e_health_model.yaml"
    expect_results = "tests/e2e_expected_outputs/test_e2e_movement_model.csv"

    delta_config = {
        "movement_model": {
            "__type__": "default_movement_model.DefaultMovementModel",
            "weighted_sampling": True,
            "home_activity": "Home_Activity",
            "location_closure_exemptions": {"House": ["Home_Activity"]},
            "facemask_activities": ["Other_Activity"],
            "age_groups": [0, 18, 65],
            "facemask_hesitancy": [0.9, 0.5, 0.2],
        }
    }

    scenario = _update_config(input_scenario, delta_config)
    _run_end_to_end_simulation(
        scenario, expect_results, tmp_path, request.node.originalname
    )


@pytest.mark.slow
def test_e2e_seasonal_effects_model(tmp_path, request):
    # input_scenario = "Scenarios/Test/test_e2e_health_seasonal_models.yaml"
    input_scenario = "Scenarios/Test/test_e2e_health_model.yaml"
    expect_results = (
        "tests/e2e_expected_outputs/test_e2e_seasonal_effects_model.csv"
    )

    delta_config = {
        "seasonal_effects_model": {
            "__type__": "default_seasonal_effects_model.DefaultSeasonalEffectsModel",
            # If the following is null, as opposed to a filepath, all multipliers are set equal to 1:
            "seasonal_multiplier_by_region_by_month": None,
            "out_of_season_multiplier": 1.0,
        }
    }

    scenario = _update_config(input_scenario, delta_config)
    _run_end_to_end_simulation(
        scenario, expect_results, tmp_path, request.node.originalname
    )


# @pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_testing_and_contact_tracing_model(tmp_path, request):
    # input_scenario = "Scenarios/Test/test_e2e_testing_and_contact_tracing_model.yaml"
    input_scenario = "Scenarios/Test/test_e2e_health_model.yaml"
    expect_results = "tests/e2e_expected_outputs/test_e2e_testing_and_contact_tracing_model.csv"

    delta_config = {
        "testing_and_contact_tracing_model": {
            "__type__": "default_testing_and_contact_tracing_model.DefaultTestingAndContactTracingModel",
            "quarantine_period_days": 14,
            "symptomatic_disease_treshold": 0.2,
            "prob_quarantine_with_symptoms_without_test": 0.0,
            "prob_quarantine_with_contact_without_test": 0.0,
            "test_threshold": 0.4,
            "test_false_negative": 0.01,
            "max_regular_contacts_to_test": 10,
        }
    }

    scenario = _update_config(input_scenario, delta_config)
    _run_end_to_end_simulation(
        scenario, expect_results, tmp_path, request.node.originalname
    )


# @pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_travel_model(tmp_path, request):
    # input_scenario = "Scenarios/Test/test_e2e_travel_model.yaml"
    input_scenario = "Scenarios/Test/test_e2e_health_model.yaml"
    expect_results = "tests/e2e_expected_outputs/test_e2e_travel_model.csv"

    delta_config = {
        "travel_model": {
            "__type__": "default_travel_model.DefaultTravelModel",
            "travel_transmission_multiplier": 1.0,
            "interpolation": 0.0,
        }
    }

    scenario = _update_config(input_scenario, delta_config)
    _run_end_to_end_simulation(
        scenario, expect_results, tmp_path, request.node.originalname
    )


@pytest.mark.slow
def test_e2e_vaccination_model(tmp_path, request):
    # input_scenario = "Scenarios/Test/test_e2e_vaccination_model.yaml"
    input_scenario = "Scenarios/Test/test_e2e_health_model.yaml"
    expect_results = "tests/e2e_expected_outputs/test_e2e_vaccination_model.csv"

    delta_config = {
        "vaccination_model": {
            "__type__": "default_vaccination_model.DefaultVaccinationModel",
            "number_of_vaccines": 2,
            "age_groups": [0, 18, 65],
            "vaccine_hesitancy": [0.8, 0.4, 0.1],
            "booster_waiting_time_days": 50,
            "vaccines": {
                "vaccine_0": {
                    "rho_immunity_failure": [
                        [[-1, 0, 14], [[1.0, 0.0], [0.5, 0.0], [1.0, 0.0]]],
                        [[-1, 0, 14], [[1.0, 0.0], [0.5, 0.0], [1.0, 0.0]]],
                    ],
                    "sigma_immunity_failure": [
                        [[-1, 0, 14], [1.0, 0.5, 1.0]],
                        [[-1, 0, 14], [1.0, 0.5, 1.0]],
                    ],
                },
                "vaccine_1": {
                    "rho_immunity_failure": [
                        [[-1, 0, 14], [[1.0, 0.0], [0.5, 0.0], [1.0, 0.0]]],
                        [[-1, 0, 14], [[1.0, 0.0], [0.5, 0.0], [1.0, 0.0]]],
                    ],
                    "sigma_immunity_failure": [
                        [[-1, 0, 14], [1.0, 0.5, 1.0]],
                        [[-1, 0, 14], [1.0, 0.5, 1.0]],
                    ],
                },
            },
        }
    }

    scenario = _update_config(input_scenario, delta_config)
    _run_end_to_end_simulation(
        scenario, expect_results, tmp_path, request.node.originalname
    )
