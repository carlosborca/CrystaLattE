"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import subprocess

def test_psi4api_ammonia():
    """Checks the psi4api mode with the ammonia crystal."""

    # Execute the main function of crystalatte and retrieve the N-mers dictionary.
    nmers, cle = crystalatte.main(
            cif_input="crystalatte/data/cif/Ammonia.cif", 
            cif_output="crystalatte/data/cif/Ammonia.xyz", 
            cif_a=3, 
            cif_b=3, 
            cif_c=3, 
            nmers_up_to=5, 
            r_cut_com=6.5, 
            r_cut_monomer=3.5, 
            r_cut_dimer=2.6, 
            r_cut_trimer=3.7, 
            r_cut_tetramer=3.7, 
            r_cut_pentamer=6.1, 
            cle_run_type=["psi4api"], 
            psi4_method="HF/STO-3G", 
            psi4_bsse="nocp", 
            psi4_memory="500 MB", 
            verbose=2
            )

    # For debugging.
    #import pprint
    #pprint.pprint(nmers)

    # Test the number of N-mers.
    number_mono  = len([k for k in nmers.keys() if k.startswith("1mer-")])
    number_di    = len([k for k in nmers.keys() if k.startswith("2mer-")])
    number_tri   = len([k for k in nmers.keys() if k.startswith("3mer-")])
    number_tetra = len([k for k in nmers.keys() if k.startswith("4mer-")])
    number_penta = len([k for k in nmers.keys() if k.startswith("5mer-")])
    
    assert compare(7, number_mono,  "Number of Monomers: ")
    assert compare(1, number_di,    "Number of Monomers: ")
    assert compare(2, number_tri,   "Number of Monomers: ")
    assert compare(1, number_tetra, "Number of Monomers: ")
    assert compare(1, number_penta, "Number of Monomers: ")

    # Test the number of atoms per monomer in each N-mer.
    for k, v in nmers.items():
        n = int(k[0])
        ref = [4]*n
        if n != 1:
            assert compare(ref, v["atoms_per_monomer"])

    # Test replicas for each N-mer.
    assert compare_values(6, nmers["2mer-0+1"]["replicas"],       "2mer-0+1 Replicas: ")
    assert compare_values(3, nmers["3mer-0+1+2"]["replicas"],     "3mer-0+1+2 Replicas: ")
    assert compare_values(3, nmers["3mer-0+1+5"]["replicas"],     "3mer-0+1+5 Replicas: ")
    assert compare_values(1, nmers["4mer-0+1+2+3"]["replicas"],   "4mer-0+1+2+3 Replicas: ")
    assert compare_values(3, nmers["5mer-0+1+2+3+4"]["replicas"], "5mer-0+1+2+3+4 Replicas: ")

    # Test the non-additive many-body energies.
    assert compare_values(0.00033410465863426,    nmers["2mer-0+1"]["nambe"],       atol=1.e-8)
    assert compare_values(1.3679718790626794e-05, nmers["3mer-0+1+2"]["nambe"],     atol=1.e-8)
    assert compare_values(-4.049166795994097e-05, nmers["3mer-0+1+5"]["nambe"],     atol=1.e-8)
    assert compare_values(2.4491332624165807e-06, nmers["4mer-0+1+2+3"]["nambe"],   atol=1.e-8)
    assert compare_values(-3.28520798120735e-07,  nmers["5mer-0+1+2+3+4"]["nambe"], atol=1.e-8)

    # Test contributions for each N-mer.
    assert compare_values(0.0010023139755617194,   nmers["2mer-0+1"]["contrib"],       atol=1.e-8)
    assert compare_values(1.3679718705361665e-05,  nmers["3mer-0+1+2"]["contrib"],     atol=1.e-8)
    assert compare_values(-4.049166741992849e-05,  nmers["3mer-0+1+5"]["contrib"],     atol=1.e-8)
    assert compare_values(6.122833582367093e-07,   nmers["4mer-0+1+2+3"]["contrib"],   atol=1.e-8)
    assert compare_values(-1.9695794435392598e-07, nmers["5mer-0+1+2+3+4"]["contrib"], atol=1.e-8)

    # Test the nuclear repulsion energies for each N-mer.
    assert compare_values(38.18102741857883,  nmers["2mer-0+1"]["nre"],       atol=1.e-8)
    assert compare_values(79.5638545934053,   nmers["3mer-0+1+2"]["nre"],     atol=1.e-8)
    assert compare_values(77.87320866497402,  nmers["3mer-0+1+5"]["nre"],     atol=1.e-8)
    assert compare_values(136.3717727214002,  nmers["4mer-0+1+2+3"]["nre"],   atol=1.e-8)
    assert compare_values(194.64346937372983, nmers["5mer-0+1+2+3+4"]["nre"], atol=1.e-8)

    # Test the COM-based priority for each N-mer.
    assert compare_values(0.0025317103172576337,  nmers["2mer-0+1"]["priority_com"],       atol=1.e-5)
    assert compare_values(2.3782189167926454e-08, nmers["3mer-0+1+2"]["priority_com"],     atol=1.e-10)
    assert compare_values(1.6227141916361213e-08, nmers["3mer-0+1+5"]["priority_com"],     atol=1.e-10)
    assert compare_values(8.289215938603657e-16,  nmers["4mer-0+1+2+3"]["priority_com"],   atol=1.e-18)
    assert compare_values(3.2435020503790996e-27, nmers["5mer-0+1+2+3+4"]["priority_com"], atol=1.e-29)

    # Test the N-mer-separation-based priority for each N-mer.
    assert compare_values(0.008552208792327147,  nmers["2mer-0+1"]["priority_min"],       atol=1.e-5)
    assert compare_values(7.570854620127571e-07, nmers["3mer-0+1+2"]["priority_min"],     atol=1.e-9)
    assert compare_values(6.255109048743461e-07, nmers["3mer-0+1+5"]["priority_min"],     atol=1.e-9)
    assert compare_values(6.937449498781642e-13, nmers["4mer-0+1+2+3"]["priority_min"],   atol=1.e-15)
    assert compare_values(1.441820074805934e-22, nmers["5mer-0+1+2+3+4"]["priority_min"], atol=1.e-24)

    # Test the N-mer-cutoff-based priority for each N-mer.
    assert compare_values(1.7489244e-03, nmers["2mer-0+1"]["priority_cutoff"],       atol=1.e-10)
    assert compare_values(6.2551090e-07, nmers["3mer-0+1+2"]["priority_cutoff"],     atol=1.e-14)
    assert compare_values(6.2551090e-07, nmers["3mer-0+1+5"]["priority_cutoff"],     atol=1.e-14)
    assert compare_values(9.3558713e-12, nmers["4mer-0+1+2+3"]["priority_cutoff"],   atol=1.e-19)
    assert compare_values(4.0052011e-25, nmers["5mer-0+1+2+3+4"]["priority_cutoff"], atol=1.e-32)

    # Test the N-mer cutoffs for each N-mer.
    assert compare_values(4.889981712498321, min(nmers["2mer-0+1"]["min_monomer_separations"]),       atol=1.e-8)
    assert compare_values(4.588498152619043, min(nmers["3mer-0+1+2"]["min_monomer_separations"]),     atol=1.e-8)
    assert compare_values(4.889981712498321, min(nmers["3mer-0+1+5"]["min_monomer_separations"]),     atol=1.e-8)
    assert compare_values(4.588498152619042, min(nmers["4mer-0+1+2+3"]["min_monomer_separations"]),   atol=1.e-8)
    assert compare_values(4.588498152619042, min(nmers["5mer-0+1+2+3+4"]["min_monomer_separations"]), atol=1.e-8)

    # Test the COM cutoffs for each N-mer.
    assert compare_values(7.337171371615357, min(nmers["2mer-0+1"]["com_monomer_separations"]),       atol=1.e-8)
    assert compare_values(6.459398295191971, min(nmers["3mer-0+1+2"]["com_monomer_separations"]),     atol=1.e-8)
    assert compare_values(7.337171371615357, min(nmers["3mer-0+1+5"]["com_monomer_separations"]),     atol=1.e-8)
    assert compare_values(6.459398295191971, min(nmers["4mer-0+1+2+3"]["com_monomer_separations"]),   atol=1.e-8)
    assert compare_values(6.459398295191971, min(nmers["5mer-0+1+2+3+4"]["com_monomer_separations"]), atol=1.e-8)

    # Test the fist eight eigenvalues of the chemical similarity matrix for each N-mer
    ref_eigenvalues = {
        "2mer-0+1" :       [61.0577465, 47.27888919,   0.47526893,  0.24647406,  0.19289154,  0.17243445, 0.15959053,  0.1341196],
        "3mer-0+1+2" :     [68.7489044, 47.64038823,  46.16025,     0.54602448,  0.27553,     0.24654741, 0.19155801,  0.18623255],
        "3mer-0+1+5" :     [67.98553018, 47.27916218, 47.27916218,  0.56261868,  0.24797375,  0.24797375, 0.21507249, 0.17507582],
        "4mer-0+1+2+3" :   [76.65031347, 47.80875755, 46.16191286, 46.16191286,  0.59866536,  0.27665228,  0.27665228,  0.25452003],
        "5mer-0+1+2+3+4" : [81.63568805, 50.70144021, 47.31855835, 46.16210417, 45.17840185,  0.64552605,  0.33151774,  0.28139148],
    }

    assert compare_values(ref_eigenvalues["2mer-0+1"], nmers["2mer-0+1"]["chsev"][:8], atol=1e-7)
    assert compare_values(ref_eigenvalues["3mer-0+1+2"], nmers["3mer-0+1+2"]["chsev"][:8], atol=1e-7)
    assert compare_values(ref_eigenvalues["3mer-0+1+5"], nmers["3mer-0+1+5"]["chsev"][:8], atol=1e-7)
    assert compare_values(ref_eigenvalues["4mer-0+1+2+3"], nmers["4mer-0+1+2+3"]["chsev"][:8], atol=1e-7)
    assert compare_values(ref_eigenvalues["5mer-0+1+2+3+4"], nmers["5mer-0+1+2+3+4"]["chsev"][:8], atol=1e-7)

    # Test the crystal lattice energy.
    assert compare_values(0.00097592, cle, atol=1.e-8)

    # Clean-up generated test files.
    subprocess.call(["rm", "-r", "crystalatte/data/cif/Ammonia"])
    subprocess.call(["rm", "crystalatte/data/cif/Ammonia.xyz"])
