"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest

def test_psz_print_header():
    """.
    """

    crystalatte.psz_print_header(2)
