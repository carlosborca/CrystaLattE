"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import time

def test_psz_print_results():
    """.
    """

    results = [] 
    results.append("2mer-0+1                   |   0.87719163 |    6 |   2.63157490 |    2.63157490 | 1.748924e-03 |  2.588")
    results.append("3mer-0+1+2                 |   0.03591606 |    3 |   0.03591606 |    2.66749096 | 6.255109e-07 |  2.428  2.588  2.588")
    results.append("3mer-0+1+5                 |  -0.10631057 |    3 |  -0.10631057 |    2.56118038 | 6.255109e-07 |  2.588  2.588  2.588")
    results.append("4mer-0+1+2+3               |   0.00643089 |    1 |   0.00160772 |    2.56278811 | 9.355871e-12 |  2.428  2.428  2.428  2.588  2.588  2.588")
    results.append("5mer-0+1+2+3+4             |  -0.00086382 |    3 |  -0.00051829 |    2.56226981 | 4.005201e-25 |  2.428  2.428  2.428  2.588  2.588  2.588  2.588  2.588  3.945  5.012")

    crystal_lattice_energy = 2.56226981

    crystalatte.psz_print_results(results, crystal_lattice_energy, 2)
