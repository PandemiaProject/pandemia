from pandemia.components.policy_maker_model.default_policy_maker_model import (
    DefaultPolicyMakerModel,
)
from pandemia.config import Config
from helper_list_models import get_clock
import pytest
from pandemia.world import region
from pandemia.random_tools import Random


@pytest.mark.xfail(
    reason="See issue https://github.com/PandemiaProject/pandemia/issues/91"
)
def test_default_policy_model_input_file(tmp_path):
    """
    Tests for the DefaultPolicyMakerModel. The model takes a directory of input files as a parameter.
    What does/should happen if the file:
    - does not exist?
    - is empty?
    - is malformed/unreadable?
    """
    # Policy file is guaranteed to not exist (because we can trust that `tmp_path` is empty)
    policy_file_path = tmp_path / "test_vector_region.csv"

    conf_dict = {
        "__type__": "default_policy_maker_model.DefaultPolicyMakerModel",
        "policy_data_filepath": str(policy_file_path),
    }

    conf = Config(_dict=conf_dict)

    dpmm = DefaultPolicyMakerModel(
        config=conf,
        scale_factor=1,
        clock=get_clock(),
        number_of_regions=4,
        number_of_vaccines=2,
        age_groups=[0, 18, 65],
    )

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

    # Check that an error is raised when the file does not exist:
    with pytest.raises(ValueError, match="does not exist"):
        dpmm.vectorize_component(vr)

    # Check that an error is when an empty file is used:
    policy_file_path.touch()
    with pytest.raises(ValueError, match="contents are not readable"):
        dpmm.vectorize_component(vr)

    # Check that an error is when an empty file is used:
    with open(policy_file_path, mode="wt") as pf:
        pf.write("test that should not parse")

    with pytest.raises(ValueError, match="contents are not readable"):
        dpmm.vectorize_component(vr)


@pytest.mark.skip("Test is not implemented")
def test_default_seasonal_multiplier_by_region_by_month():
    """
    Test the inputs to the DefaultSeasonalEffectsModel
    - what does/should happen if `seasonal_multiplier_by_region_by_month is None`?
      (which appear to be a valid condition)
    - what does/should happen if `seasonal_multiplier_by_region_by_month` points to an non-existant file?
    - what does/should happen if there is a mismatch between regions in file and regions defined elsewhere
    """
    pytest.fail()


@pytest.mark.skip("Test is not implemented")
def test_default_vaccination_model_inputs():
    """
    Test the inputs to the DefaultVaccinationModel.
    In would be instructive to have tests around the values for `number_of_strains` and `number_of_rho_immunity_outcomes`.

    Why are these required given the values can be calculated from values in the health model? Is this
    to ensure that the user input is self-consistent? If so great, but it would good to have more user-friendly
    error messages if so.
    (In general is it better to use ValueError for things that are errors in user input, and
    use assert/AssertionError for unittests and self-consistency of the software - not the input.
    """
    pytest.fail()
