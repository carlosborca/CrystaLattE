# <img align="center" src="https://github.com/carlosborca/CrystaLattE/blob/master/media/logo/Logo.png" height=260>

Automated calculation of crystal lattice energies with the many-body expansion.

| Category | Badges |
|-------------|-------------|
| **Status** | [![Travis Build Status](https://travis-ci.com/carlosborca/CrystaLattE.svg?branch=master)](https://travis-ci.org/carlosborca/CrystaLattE) [![codecov](https://codecov.io/gh/carlosborca/CrystaLattE/branch/master/graph/badge.svg)](https://codecov.io/gh/carlosborca/CrystaLattE/branch/master) [![Total alerts](https://img.shields.io/lgtm/alerts/g/carlosborca/CrystaLattE.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/carlosborca/CrystaLattE/alerts/) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/carlosborca/CrystaLattE.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/carlosborca/CrystaLattE/context:python) |
| **Foundation** | [![License](https://img.shields.io/github/license/carlosborca/CrystaLattE.svg)](https://opensource.org/licenses/LGPL-3.0) [![GitHub Top Languages](https://img.shields.io/github/languages/top/carlosborca/CrystaLattE)](https://github.com/carlosborca/CrystaLattE/) |
| **GitHub Info** | [![GitHub Code Size](https://img.shields.io/github/languages/code-size/carlosborca/CrystaLattE)](https://github.com/carlosborca/CrystaLattE/) [![GitHub Commits per Month](https://img.shields.io/github/commit-activity/m/carlosborca/CrystaLattE)](https://github.com/carlosborca/CrystaLattE/) [![GitHub Last Commit](https://img.shields.io/github/last-commit/carlosborca/CrystaLattE)](https://github.com/carlosborca/CrystaLattE/) |

## Overview

This is the overview.

## General Usage information

Some instructions here.

### Installation

Minimal set of commands to install CrystaLattE on Linux, MacOS, or Windows (with the Windows Subsystem for Linux). Last tested on 4 October 2019:

#### Install miniconda:

Get the installer from the Anaconda website:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Check the installer (Optional):

```
sha256sum Miniconda3-latest-Linux-x86_64.sh
```

Run the installer following the on-screen instructions:

```
bash Miniconda3-latest-Linux-x86_64.sh
```

After the installation is complete, restart the shell.

Disable automatic activation of the _base_ conda environment _(Optional)_:

```
conda config --set auto_activate_base false
```

Test the conda installation:

```
conda list
```

#### Create a conda environment for CrystaLattE:

```
conda create -n CrystaLattE python=3.7 psi4 psi4-rt pycifrw -c psi4/label/dev -c psi4 -c conda-forge
```

Note: The `v2rdm_casscf` plugin of PSI4 can cause errors when trying to import PSI4 from CrystaLattE. In that case, remove it.

```
conda remove v2rdm_casscf
```

#### Activate the CrystaLattE conda environment:

```
conda activate CrystaLattE
```

#### Clone CrystaLattE from its GitHub repository:

```
git clone https://github.com/carlosborca/CrystaLattE.git
```

#### Go to the root directory and run the test suite:

```
cd CrystaLattE/; pytest
```

### How to run CrystaLattE

#### Preparing the options input

Need a CIF file and an Options input file for CrystaLattE.

```
# CrystaLattE Input Template:

# Blank lines and lines starting with the hash character are omitted.
# Padding blank spaces are also ignored.

# The format is:
# Keywords      = Value

cif_input       = Benzene.cif
cif_output      = Benzene.xyz
cif_a           = 5
cif_b           = 5
cif_c           = 5
nmers_up_to     = 5
r_cut_com       = 8.4
r_cut_monomer   = 4.43
r_cut_dimer     = 2.6
r_cut_trimer    = 2.7
r_cut_tetramer  = 2.7
r_cut_pentamer  = 5.1
cle_run_type    = psi4api + quiet
psi4_method     = MP2/aug-cc-pV[TQ]Z + D:FNO-CCDS(T)/aug-cc-pVDZ
psi4_bsse       = cp
psi4_memory     = 500 MB
verbose         = 2
```

#### Running CrystaLattE

```
crystalatte.py input.cle
```

#### Copyright

Copyright (c) 2019, Carlos H. Borca


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.0.
