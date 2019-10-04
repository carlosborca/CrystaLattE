# <img src="https://github.com/carlosborca/CrystaLattE/blob/master/media/logo/Logo.png" height=260>

| **Status** | [![Travis Build Status](https://travis-ci.com/carlosborca/CrystaLattE.svg?branch=master)](https://travis-ci.org/carlosborca/CrystaLattE) [![codecov](https://codecov.io/gh/carlosborca/CrystaLattE/branch/master/graph/badge.svg)](https://codecov.io/gh/carlosborca/CrystaLattE/branch/master) |
|-------------|:-------------:|
| **Foundation** | [![license](https://img.shields.io/github/license/psi4/psi4.svg)](https://opensource.org/licenses/LGPL-3.0) [![platforms](https://img.shields.io/badge/Platforms-Linux%2C%20MacOS%2C%20Windows%20WSL-orange.svg)](http://psicode.org/psi4manual/master/introduction.html#supported-systems) [![python](https://img.shields.io/badge/python-3.6+-blue.svg)](http://psicode.org/psi4manual/master/introduction.html#supported-systems) |

Automated calculation of crystal lattice energies with the many-body expansion.

## Overview

This is the overview.

## Installation

Minimal set of commands to install and run CrystaLattE on Linux, MacOS, or Windows (with the Windows Subsystem for Linux). Last tested on 4 October 2019:

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


#### Copyright

Copyright (c) 2019, Carlos H. Borca


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.0.
