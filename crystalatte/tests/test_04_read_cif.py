"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest

def test_cle_read_cif():
    """Checks that the read CIF function can take a CIF and produce a
    correct dictionary containing the unit cell structural data for an
    ammonia crystal."""

    ref = {
        '_cell_length_a': 5.1305,
        '_cell_length_b': 5.1305,
        '_cell_length_c': 5.1305,
        '_cell_angle_alpha': 90.0,
        '_cell_angle_beta': 90.0,
        '_cell_angle_gamma': 90.0,
        '_cell_volume': 135.05,
        '_symmetry_equiv_pos_as_xyz': ['x,y,z'],
        '_atom_site_label': ['N', 'N', 'N', 'N', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],
        '_atom_site_fract_x': [0.2107, 0.7107, 0.7893, 0.2893, 0.3689, 0.2671, 0.1159, 0.8689, 0.7671, 0.6159, 0.7329, 0.8841, 0.3841, 0.1311, 0.2329, 0.6311],
        '_atom_site_fract_y': [0.2107, 0.2893, 0.7107, 0.7893, 0.2671, 0.1159, 0.3689, 0.2329, 0.3841, 0.1311, 0.6159, 0.8689, 0.6311, 0.7329, 0.8841, 0.7671],
        '_atom_site_fract_z': [0.2107, 0.7893, 0.2893, 0.7107, 0.1159, 0.3689, 0.2671, 0.8841, 0.6311, 0.7329, 0.1311, 0.2329, 0.7671, 0.6159, 0.8689, 0.3841]
        }

    data = crystalatte.read_cif("crystalatte/data/cif/Ammonia.cif")

    assert compare(ref, data)
