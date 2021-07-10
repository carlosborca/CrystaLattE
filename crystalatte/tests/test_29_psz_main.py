"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import os
import pytest
import subprocess

def test_psz_main():
    """Checks that the main function on the psithonyzer script returns
    the correct values for precomputed psithon outputs."""
    
    kwargs = {'com_mode': False, 'sort_by_avg_com_dist': False, 'sort_by_nmer_cutoff': False, 'src_directory': "crystalatte/data/out"}
    results, crystal_lattice_energy = crystalatte.psz_main(2, **kwargs)


    a = []
    a.append("2mer-0+1                   |   0.87719178 |    6 |   2.63157535 |    2.63157535 | 1.748924e-03 |  2.588 ")
    a.append("3mer-0+1+2                 |   0.03591565 |    3 |   0.03591565 |    2.66749100 | 6.255109e-07 |  2.428  2.588  2.588 ")
    a.append("3mer-0+1+5                 |  -0.10631085 |    3 |  -0.10631085 |    2.56118015 | 6.255109e-07 |  2.588  2.588  2.588 ")
    a.append("4mer-0+1+2+3               |   0.00643113 |    1 |   0.00160778 |    2.56278793 | 9.355871e-12 |  2.428  2.428  2.428  2.588  2.588  2.588 ")
    a.append("5mer-0+1+2+3+4             |  -0.00086359 |    3 |  -0.00051816 |    2.56226977 | 4.005201e-25 |  2.428  2.428  2.428  2.588  2.588  2.588  2.588  2.588  3.945  5.012 ")

    assert compare(a[0], results[0])
    assert compare(a[1], results[1])
    assert compare(a[2], results[2])
    assert compare(a[3], results[3])
    assert compare(a[4], results[4])
    assert compare_values(2.56226977, crystal_lattice_energy, atol=1.e-9)
    
    # Clean-up generated test files.
    subprocess.call(["rm", "Ammonia.csv"])
