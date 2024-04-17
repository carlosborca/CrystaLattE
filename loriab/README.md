Several starting points copied from Philip

* CLE repo is a clone from Philip that includes his and Zach's (and probably others') corrections: https://github.com/philipmnel/CrystaLattE.git, branch `philip`, commit `d70da56`
* `cif/acetic_acid.cif` copied from `/theoryfs2/ds/pmnelson/chem/crystals/cif/acetic_acid.cif`
* `acetic_acid_mp2_noaug.cle` copied from `/theoryfs2/ds/pmnelson/chem/crystals/acetic_acid_mp2_noaug.cle`

Directions
* Create a conda env with psi4 (brings networkx) and pycifrw. c-f can provide both.
* From repo-top-level dir, run `pytest crystalatte/` to test.
* Create any directories in the `cif_output` option before running CLE.
* The `cle_run_type    = qcmanybody` writes mol json to `cif_output` and runs run_qcengine
* Run CLE as `../crystalatte/crystalatte.py demo.cle`


