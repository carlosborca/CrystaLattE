"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
import crystalatte
import pytest

def test_noninteger_nmers_up_to_input():
    """Checks that the program prints an error message when an invalid
    output with a non-integer nmers_up_to keyword is provided."""

    with pytest.raises(SystemExit):
        with pytest.raises(IndexError):
            crystalatte.input_parser("crystalatte/data/cle/int_keywords.cle")
