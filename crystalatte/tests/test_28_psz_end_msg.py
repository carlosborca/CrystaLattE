"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import time

def test_psz_print_end_msg():
    """.
    """

    start = time.time()
    crystalatte.psz_print_end_msg(start, 2)
