"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
import crystalatte
import pytest

def test_nonfloat_r_cut_com_input():
    """Checks that the program prints an error message when an invalid
    output with a non-float 'r_cut_com' keyword is provided."""

    with pytest.raises(SystemExit):
        with pytest.raises(IndexError):
            crystalatte.input_parser("crystalatte/data/cle/float_keywords.cle")
