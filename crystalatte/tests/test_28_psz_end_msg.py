"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import time

def test_psz_print_header():
    """Checks that the psithonyzer success checker function retruns a
    True boolean value when it finds a Psi4's successful execution
    string."""

    start = time.time()
    crystalatte.psz_print_end_msg(start, 2)
