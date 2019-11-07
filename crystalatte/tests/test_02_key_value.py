"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
import crystalatte
import pytest

def test_wrong_input():
    """Test that the program exits when a wrong CLE input is provided."""

    with pytest.raises(SystemExit):
        with pytest.raises(IndexError):
            crystalatte.input_parser("crystalatte/data/cle/key_value.cle")
