"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import subprocess

@pytest.mark.parametrize("run_type,ref_cle",
                         [
                             pytest.param(["test", "timings"], 0.0, id="timings"),
                             pytest.param(["psithon"], 0.0, id="psithon"),
                             pytest.param(["psi4api"], 0.0009759168, id="psi4api")
                         ])
def test_timings_ammonia(run_type, ref_cle):
    """Checks the ammonia crystal results in various run modes."""

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
            cle_run_type=run_type, 
            psi4_method="HF/STO-3G", 
            psi4_bsse="nocp", 
            psi4_memory="500 MB", 
            verbose=2
            )

    # For debugging.
    import pprint
    pprint.pprint(nmers)

    # Test the number of N-mers.
    number_mono  = len([k for k in nmers.keys() if k.startswith("1mer-")])
    number_di    = len([k for k in nmers.keys() if k.startswith("2mer-")])
    number_tri   = len([k for k in nmers.keys() if k.startswith("3mer-")])
    number_tetra = len([k for k in nmers.keys() if k.startswith("4mer-")])
    number_penta = len([k for k in nmers.keys() if k.startswith("5mer-")])
    
    assert compare(7, number_mono,  "Number of Monomers: ")
    assert compare(1, number_di,    "Number of Dimers: ")
    assert compare(2, number_tri,   "Number of Trimers: ")
    assert compare(1, number_tetra, "Number of Tetramers: ")
    assert compare(1, number_penta, "Number of Pentamers: ")

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

    # Test the nuclear repulsion energies for each N-mer.
    assert compare_values(38.181027098543616,  nmers["2mer-0+1"]["nre"],       atol=1.e-8)
    assert compare_values(79.56385358458158,   nmers["3mer-0+1+2"]["nre"],     atol=1.e-8)
    assert compare_values(77.87321016752901,  nmers["3mer-0+1+5"]["nre"],     atol=1.e-8)
    assert compare_values(136.37176983414787,  nmers["4mer-0+1+2+3"]["nre"],   atol=1.e-8)
    assert compare_values(194.64346863543244, nmers["5mer-0+1+2+3+4"]["nre"], atol=1.e-8)

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
    assert compare_values(1.7489237787727267e-03, nmers["2mer-0+1"]["priority_cutoff"],       atol=1.e-10)
    assert compare_values(6.255103294673478e-07, nmers["3mer-0+1+2"]["priority_cutoff"],     atol=1.e-14)
    assert compare_values(6.255103294189915e-07, nmers["3mer-0+1+5"]["priority_cutoff"],     atol=1.e-14)
    assert compare_values(9.35585604123943e-12, nmers["4mer-0+1+2+3"]["priority_cutoff"],   atol=1.e-19)
    assert compare_values(4.005207215747089e-25, nmers["5mer-0+1+2+3+4"]["priority_cutoff"], atol=1.e-32)

    # Test the N-mer cutoffs for each N-mer.
    assert compare_values(4.889982212459292, min(nmers["2mer-0+1"]["min_monomer_separations"]),       atol=1.e-8)
    assert compare_values(4.588498416229297, min(nmers["3mer-0+1+2"]["min_monomer_separations"]),     atol=1.e-8)
    assert compare_values(4.889982212459292, min(nmers["3mer-0+1+5"]["min_monomer_separations"]),     atol=1.e-8)
    assert compare_values(4.588498416229297, min(nmers["4mer-0+1+2+3"]["min_monomer_separations"]),   atol=1.e-8)
    assert compare_values(4.588498416229297, min(nmers["5mer-0+1+2+3+4"]["min_monomer_separations"]), atol=1.e-8)

    # Test the COM cutoffs for each N-mer.
    assert compare_values(7.337170476925919, min(nmers["2mer-0+1"]["com_monomer_separations"]),       atol=1.e-8)
    assert compare_values(6.459398895035611, min(nmers["3mer-0+1+2"]["com_monomer_separations"]),     atol=1.e-8)
    assert compare_values(7.337170476925917, min(nmers["3mer-0+1+5"]["com_monomer_separations"]),     atol=1.e-8)
    assert compare_values(6.459398895035611, min(nmers["4mer-0+1+2+3"]["com_monomer_separations"]),   atol=1.e-8)
    assert compare_values(6.459398895035611, min(nmers["5mer-0+1+2+3+4"]["com_monomer_separations"]), atol=1.e-8)

    # Test the fist eight eigenvalues of the chemical similarity matrix for each N-mer
    ref_eigenvalues = {
        "2mer-0+1" :       [61.057747426825, 47.278888122977, 0.475268780429, 0.246473958505, 0.192891650536, 0.172434560462, 0.159590615089, 0.134119684833],
        "3mer-0+1+2" :     [68.74890504188, 47.640386591887, 46.160250779761, 0.546024355783, 0.275529905072, 0.246547337798, 0.191558119533, 0.186232649809],
        "3mer-0+1+5" :     [67.985532101019, 47.279161107153, 47.279161107153, 0.562618519811, 0.247973654359, 0.247973654359, 0.215072609408, 0.17507592859],
        "4mer-0+1+2+3" :   [76.65031350828, 47.808755663003, 46.16191363632, 46.16191363632, 0.598665256802, 0.276652190605, 0.276652190605, 0.254519979836],
        "5mer-0+1+2+3+4" : [81.635688941068, 50.701439488015, 47.318557217411, 46.162104945763, 45.178401675932, 0.645525934722, 0.331517648871, 0.281391408722],
    }

    assert compare_values(ref_eigenvalues["2mer-0+1"], nmers["2mer-0+1"]["chsev"][:8], atol=1e-7)
    assert compare_values(ref_eigenvalues["3mer-0+1+2"], nmers["3mer-0+1+2"]["chsev"][:8], atol=1e-7)
    assert compare_values(ref_eigenvalues["3mer-0+1+5"], nmers["3mer-0+1+5"]["chsev"][:8], atol=1e-7)
    assert compare_values(ref_eigenvalues["4mer-0+1+2+3"], nmers["4mer-0+1+2+3"]["chsev"][:8], atol=1e-7)
    assert compare_values(ref_eigenvalues["5mer-0+1+2+3+4"], nmers["5mer-0+1+2+3+4"]["chsev"][:8], atol=1e-7)

    # Test the crystal lattice energy.
    assert compare_values(ref_cle, cle, atol=1.e-8)
    
    # Clean-up generated test files.
    subprocess.call(["rm", "-r", "crystalatte/data/cif/Ammonia"])
    files_to_clean = ["crystalatte/data/cif/Ammonia.xyz",
                      "2mer-0+1.dat",
                      "3mer-0+1+2.dat",
                      "3mer-0+1+5.dat",
                      "4mer-0+1+2+3.dat",
                      "5mer-0+1+2+3+4.dat"]
    for f in files_to_clean:
        subprocess.call(["rm", f])
