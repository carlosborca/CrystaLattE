"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import subprocess

@pytest.mark.parametrize("psi4_version", ["old", "new"])
def test_psz_psz_get_nmer_data(psi4_version):
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
        if psi4_version == "old":
            f1.write("      n-Body     Total Energy [Eh]       I.E. [kcal/mol]      Delta [kcal/mol]\n")
            f1.write("           1     -166.349896077773        0.000000000000        0.000000000000\n")
            f1.write("           2     -166.352512228082       -1.641659103685       -1.641659103685\n")
            f1.write("           3     -166.352498548533       -1.633075056645        0.008584047040\n")
        elif psi4_version == "new":
            f1.write("      n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy\n")
            f1.write("                 [Eh]                    [Eh]                  [kcal/mol]            [Eh]                  [kcal/mol]\n")
            f1.write("           1  N/A                         0.000000000000        0.000000000000        0.000000000000        0.000000000000\n")
            f1.write("           2  N/A                         0.000000000000       -1.641659103685        0.000000000000       -1.641659103685\n")
            f1.write("FULL/RTN   3  N/A                         0.000000000000       -1.633075056645        0.000000000000        0.008584047040\n")
        f1.write("\n")
        f1.write("Psi4 exiting successfully. Buy a developer a beer!")

    ret = crystalatte.psz_get_nmer_data("3mer-0+1+2.out", 2)

    assert compare("Ammonia", ret["sc_xyz"])
    assert compare("3mer-0+1+2", ret["key"])
    assert compare(3, ret["number_of_monomers"])
    assert compare(3, ret["replicas"])
    assert compare(6.255109050478e-07, ret["priority_cutoff"])
    assert compare(7.570854620128e-07, ret["priority_min"])
    assert compare(2.378218916793e-08, ret["priority_com"])
    assert compare(['2.428', '2.588', '2.588'], ret["min_monomer_separations"])
    assert compare(['3.418', '3.883', '3.883'], ret["com_monomer_separations"])
    assert compare(79.5638545934053, ret["nre"])
    assert compare(0.03591565281536, ret["nambe"])

    # Clean-up generated test files.
    subprocess.call(["rm", "3mer-0+1+2.out"])
    subprocess.call(["rm", "3mer-0+1+5.out"])

@pytest.mark.skip(reason="Legacy output format no longer supported.")
@pytest.mark.parametrize("psi4_version", ["old", "new"])
def test_psz_psz_get_nmer_data_legacy(psi4_version):
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
        if psi4_version == "old":
            f2.write("   n-Body     Total Energy [Eh]       I.E. [kcal/mol]      Delta [kcal/mol]\n")
            f2.write("        1     -166.349896077479        0.000000000000        0.000000000000\n")
            f2.write("        2     -166.348893763567        0.628961475807        0.628961475807\n")
            f2.write("        3     -166.348934255125        0.603552639217       -0.025408836590\n")
        elif psi4_version == "new":
            f2.write("      n-Body     Total Energy            Interaction Energy                          N-body Contribution to Interaction Energy\n")
            f2.write("                 [Eh]                    [Eh]                  [kcal/mol]            [Eh]                  [kcal/mol]\n")
            f2.write("           1  N/A                         0.000000000000        0.000000000000        0.000000000000        0.000000000000\n")
            f2.write("           2  N/A                         0.000000000000        0.628961475807        0.000000000000        0.628961475807\n")
            f2.write("FULL/RTN   3  N/A                         0.000000000000        0.603552639217        0.000000000000       -0.025408836590\n")
        f2.write("\n")
        f2.write("Psi4 exiting successfully. Buy a developer a beer!")

    ret = crystalatte.psz_get_nmer_data("3mer-0+1+5.out", 2)

    assert compare(None, ret["sc_xyz"])
    assert compare("3mer-0+1+5", ret["key"])
    assert compare(3, ret["number_of_monomers"])
    assert compare(3, ret["replicas"])
    assert compare(6.25510905e-07, ret["priority_min"])
    assert compare(1.62271419e-08, ret["priority_com"])
    assert compare(['2.588', '2.588', '2.588'], ret["min_monomer_separations"])
    assert compare(['3.883', '3.883', '3.883'], ret["com_monomer_separations"])
    assert compare(None, ret["nre"])
    assert compare(-0.10631057229256001, ret["nambe"])
    assert compare_values(1.920089851168e-04, ret["p_cutoff"],  atol=1.e-15)

    # Clean-up generated test files.
    subprocess.call(["rm", "3mer-0+1+2.out"])
    subprocess.call(["rm", "3mer-0+1+5.out"])
