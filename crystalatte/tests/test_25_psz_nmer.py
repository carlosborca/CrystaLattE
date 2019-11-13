"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import subprocess

def test_psz_psz_get_nmer_data():
    """Checks that the returns of the psithonyzer function to retrieve
    information from the N-mer output work correctly."""

    # Newer psithon format test.
    with open("3mer-0+1+2.out", "w") as f1:
        f1.write("# This is a dummy Psi4 output test file.\n")
        f1.write("# PSI4 file produced by CrystaLattE\n")
        f1.write("# Generated from:               Some/Path/To/Ammonia\n")
        f1.write("# Psithon input for N-mer:      3mer-0+1+2\n")
        f1.write("# Number of atoms per monomer:  [4, 4, 4]\n")
        f1.write("# Number of replicas:           3\n")
        f1.write("# COM priority:                 2.378218916793e-08\n")
        f1.write("# Minimum COM separations:      3.418  3.883  3.883\n")
        f1.write("# Separation priority:          7.570854620128e-07\n")
        f1.write("# Minimum monomer separations:  2.428  2.588  2.588\n")
        f1.write("# Cutoff priority:              6.255109050478e-07\n")
        f1.write("# Nuclear repulsion energy:     79.5638545934053\n")
        f1.write("\n")
        f1.write("   ==> N-Body: Non-Counterpoise Corrected (NoCP) energies <==\n")
        f1.write("\n")
        f1.write("      n-Body     Total Energy [Eh]       I.E. [kcal/mol]      Delta [kcal/mol]\n")
        f1.write("           1     -166.349896077773        0.000000000000        0.000000000000\n")
        f1.write("           2     -166.352512228082       -1.641659103685       -1.641659103685\n")
        f1.write("           3     -166.352498548533       -1.633075056645        0.008584047040\n")
        f1.write("\n")
        f1.write("Psi4 exiting successfully. Buy a developer a beer!")

    sc_xyz, key, number_of_monomers, replicas, p_cutoff, p_min, p_com, min_mon_seps, com_mon_seps, nre, n_body_energy = crystalatte.psz_get_nmer_data("3mer-0+1+2.out", 2)

    assert compare("Ammonia", sc_xyz)
    assert compare("3mer-0+1+2", key)
    assert compare(3, number_of_monomers)
    assert compare(3, replicas)
    assert compare(6.255109050478e-07, p_cutoff)
    assert compare(7.570854620128e-07, p_min)
    assert compare(2.378218916793e-08, p_com)
    assert compare(['2.428', '2.588', '2.588'], min_mon_seps)
    assert compare(['3.418', '3.883', '3.883'], com_mon_seps)
    assert compare(79.5638545934053, nre)
    assert compare(0.03591565281536, n_body_energy)

    # Legacy psithon format test.
    with open("3mer-0+1+5.out", "w") as f2:

        f2.write("# PSI4 file produced by CrystaLattE\n")
        f2.write("# Psithon input for N-mer:      3mer-0+1+5\n")
        f2.write("# Number of replicas:           3\n")
        f2.write("# Priority index for input:     6.25510905e-07\n")
        f2.write("# Minimum monomer separations:  2.588  2.588  2.588 \n")
        f2.write("# COM priority index for input: 1.62271419e-08\n")
        f2.write("# Minimum COM separations:      3.883  3.883  3.883 \n")
        f2.write("\n")
        f2.write("   ==> N-Body: Non-Counterpoise Corrected (NoCP) energies <==\n")
        f2.write("\n")
        f2.write("   n-Body     Total Energy [Eh]       I.E. [kcal/mol]      Delta [kcal/mol]\n")
        f2.write("        1     -166.349896077479        0.000000000000        0.000000000000\n")
        f2.write("        2     -166.348893763567        0.628961475807        0.628961475807\n")
        f2.write("        3     -166.348934255125        0.603552639217       -0.025408836590\n")
        f2.write("\n")
        f2.write("Psi4 exiting successfully. Buy a developer a beer!")

    sc_xyz, key, number_of_monomers, replicas, p_cutoff, p_min, p_com, min_mon_seps, com_mon_seps, nre, n_body_energy = crystalatte.psz_get_nmer_data("3mer-0+1+5.out", 2)

    assert compare(None, sc_xyz)
    assert compare("3mer-0+1+5", key)
    assert compare(3, number_of_monomers)
    assert compare(3, replicas)
    assert compare(6.25510905e-07, p_min)
    assert compare(1.62271419e-08, p_com)
    assert compare(['2.588', '2.588', '2.588'], min_mon_seps)
    assert compare(['3.883', '3.883', '3.883'], com_mon_seps)
    assert compare(None, nre)
    assert compare(-0.10631057229256001, n_body_energy)
    assert compare_values(1.920089851168e-04, p_cutoff,  atol=1.e-15)

    # Clean-up generated test files.
    subprocess.call(["rm", "3mer-0+1+2.out"])
    subprocess.call(["rm", "3mer-0+1+5.out"])
