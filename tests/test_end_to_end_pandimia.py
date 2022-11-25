import sys
import pandas as pd
import pandemia
import pytest
import yaml

from icecream import ic

@pytest.mark.slow
def test_end_to_end(tmp_path):
    fname_input_scenario = "Scenarios/Test/test_global_config.yaml"
    fname_altered_scenario = str(tmp_path / "test_global_config.yaml")

    fname_expect_results = "tests/e2e_expected_outputs/strain_counts.csv"
    fname_actual_results = str(tmp_path / "strain_counts.csv")

    ic(fname_altered_scenario)
    ic(fname_actual_results)

    # Update the config with temporary dir for the output files
    with open(fname_input_scenario, "r") as scen:
        scenario = yaml.safe_load(scen)

    scenario["reporters"]["csv.StrainCounts"]["filename"] = fname_actual_results

    with open(fname_altered_scenario, "w") as scen:
        yaml.dump(scenario, scen)

    # Fake the commandline args and call Pandemia.main() 
    sys.argv = [pandemia.__file__, fname_altered_scenario]
    pandemia.main()

    # Actual results
    actual_df = pd.read_csv(fname_actual_results)

    # Expected results
    expected_df = pd.read_csv(fname_expect_results)

    # Check for differences in the results
    diff_df = expected_df.compare(actual_df)
    
    assert len(diff_df) == 0
    assert len(diff_df.index) == 0
    assert len(diff_df.columns) == 0
