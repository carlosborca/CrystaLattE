"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import numpy as np
import pytest
import subprocess

def test_center_supercell():
    """Checks the routine to read the supercell XYZ file and create
    Numpy arrays from it containing the center of the supercell, its
    coordinates and elements."""

    crystalatte.cif_main(fNameIn="crystalatte/data/cif/Ammonia.cif",
                         fNameOut="sc.xyz",
                         Na=1,
                         Nb=1,
                         Nc=1,
                         monomer_cutoff=0.0,
                         nmer_cutoff=0.0,
                         make_rect_box=True)

    # Execute the main function of crystalatte and retrieve the N-mers dictionary.
    scell_geom_max_coords, scell_geom, scell_elem = crystalatte.center_supercell("sc.xyz", 2)
    
    scgmc = np.array([4.84761994, 4.84761994, 4.84761994])
    sce = np.array(['N', 'N', 'N', 'N', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'])
    
    l_scg = [[-0.91910799,  3.92851196, 3.16646478],
        [ 3.92851196,  3.16646478, -0.91910799],
        [ 3.16646478, -0.91910799,  3.92851196], 
        [-1.68115516, -1.68115516, -1.68115516],
        [ 2.3947252 ,  3.71327782,  0.        ],
        [-1.46592103,  4.84761994,  4.70025154],
        [-2.45289475,  3.38169892,  2.24735679],
        [ 0.        ,  2.3947252 ,  3.71327782],
        [ 4.84761994,  4.70025154, -1.46592103],
        [ 3.38169892,  2.24735679, -2.45289475],
        [ 2.24735679, -2.45289475,  3.38169892],
        [ 3.71327782,  0.        ,  2.3947252 ],
        [ 4.70025154, -1.46592103,  4.84761994],
        [-2.60026315, -0.1473684 , -1.13434212],
        [-1.13434212, -2.60026315, -0.1473684 ],
        [-0.1473684 , -1.13434212, -2.60026315]]
    
    scg = np.array(l_scg)

    assert compare_values(scgmc, scell_geom_max_coords)
    assert compare_values(scg, scell_geom)
    assert compare(sce, scell_elem)

    # Clean-up generated test files
    subprocess.call(["rm", "sc.xyz"])
