import sys
import pandas as pd
import pandemia
import pytest
import yaml

from icecream import ic

def run_end_to_end_simulation(input_scenario_config, expected_results_csv, tmp_path):
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
    path, or relative to the the root of the repo. The test     # input_scenario = "Scenarios/Test/test_config.yaml"
    # expect_results = "tests/e2e_expected_outputs/strain_counts.csv"

    # run_end_to_end_simulation(input_scenario, expect_results, tmp_path)
es or fails based on the comparision of
    this csv file and the output of the simulation will be with this csv, using the
    `pandas.DataFrame.compare` function.
    @tmp_path A `pathlib.Path` object, representing a temporary directory to use to write outputs.
    Typically this the `tmp_path` fixture     # input_scenario = "Scenarios/Test/test_config.yaml"
    # expect_results = "tests/e2e_expected_outputs/strain_counts.csv"

    # run_end_to_end_simulation(input_scenario, expect_results, tmp_path)
ed from a calling `test_` wrapper function.
    https://docs.pytest.org/en/latest/how-to/tmp_path.html#the-tmp-path-fixture

    """
    # fname_input_scenario = "Scenarios/Test/test_global_config.yaml"
    fname_altered_scenario = str(tmp_path / "test_altered_scenario.yaml")

    # fname_expect_results = "tests/e2e_expected_outputs/strain_counts.csv"
    fname_actual_results = str(tmp_path / "actual_test_strain_counts.csv")

    ic(fname_altered_scenario)
    ic(fname_actual_results)

    # Update the config with temporary dir for the output files
    with open(input_scenario_config, "r") as scen:
        scenario = yaml.safe_load(scen)

    scenario["reporters"]["csv.StrainCounts"]["filename"] = fname_actual_results

    with open(fname_altered_scenario, "w") as scen:
        yaml.dump(scenario, scen)

    # Read the expected results
    # Do this before the main simulation, so that the test fails quickly in the case of a 
    # missing/unreadable `expected_results_csv` file.
    expected_df = pd.read_csv(expected_results_csv)

    # Fake the commandline args and call Pandemia.main() 
    sys.argv = [pandemia.__file__, fname_altered_scenario]
    pandemia.main()

    # Actual results
    actual_df = pd.read_csv(fname_actual_results)

    # Check for differences in the results
    diff_df = expected_df.compare(actual_df)
    
    assert len(diff_df) == 0
    assert len(diff_df.index) == 0
    assert len(diff_df.columns) == 0


@pytest.mark.slow
def test_end_to_end_global(tmp_path):
    input_scenario = "Scenarios/Test/test_global_config.yaml"
    # expect_results = "tests/e2e_expected_outputs/strain_counts.csv"
    expect_results = "tests/e2e_expected_outputs/strain_counts_test_global_config.csv"

    run_end_to_end_simulation(input_scenario, expect_results, tmp_path)

@pytest.mark.slow
def test_end_to_end_all_components(tmp_path):
    input_scenario = "Scenarios/Test/test_all_components.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts_test_all_components.csv"

    run_end_to_end_simulation(input_scenario, expect_results, tmp_path)

@pytest.mark.slow
def test_end_to_end_all_void(tmp_path):
    input_scenario = "Scenarios/Test/test_void_all.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts_test_void_all.csv"

    run_end_to_end_simulation(input_scenario, expect_results, tmp_path)

@pytest.mark.slow
def test_e2e_health_model(tmp_path):
    input_scenario = "Scenarios/Test/test_e2e_health_model.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts_test_e2e_health_model.csv"

    run_end_to_end_simulation(input_scenario, expect_results, tmp_path)

@pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_hospitalization_and_death_model(tmp_path):
    # input_scenario = "Scenarios/Test/test_e2e_hospitalization_and_death_model.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts.csv"

    # run_end_to_end_simulation(input_scenario, expect_results, tmp_path)

@pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_input_model(tmp_path):
    # input_scenario = "Scenarios/Test/test_e2e_input_model.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts.csv"

    # run_end_to_end_simulation(input_scenario, expect_results, tmp_path)

@pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_movement_model(tmp_path):
    # input_scenario = "Scenarios/Test/test_e2e_movement_model.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts.csv"

    # run_end_to_end_simulation(input_scenario, expect_results, tmp_path)

@pytest.mark.slow
def test_e2e_seasonal_effects_model(tmp_path):
    input_scenario = "Scenarios/Test/test_e2e_health_seasonal_models.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts_test_e2e_health_seasonal_model.csv"

    run_end_to_end_simulation(input_scenario, expect_results, tmp_path)

@pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_testing_and_contact_tracing_model(tmp_path):
    # input_scenario = "Scenarios/Test/test_e2e_testing_and_contact_tracing_model.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts.csv"

    # run_end_to_end_simulation(input_scenario, expect_results, tmp_path)

@pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_travel_model(tmp_path):
    # input_scenario = "Scenarios/Test/test_e2e_travel_model.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts.csv"

    # run_end_to_end_simulation(input_scenario, expect_results, tmp_path)

@pytest.mark.skip("Not implemented yet")
@pytest.mark.slow
def test_e2e_vaccination_model(tmp_path):
    # input_scenario = "Scenarios/Test/test_e2e_vaccination_model.yaml"
    expect_results = "tests/e2e_expected_outputs/strain_counts.csv"

    # run_end_to_end_simulation(input_scenario, expect_results, tmp_path)
