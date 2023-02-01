import pytest

@pytest.mark.skip("not required")
def test_always_pass():
    assert True

@pytest.mark.skip("not required")
def test_always_fail():
    assert False
