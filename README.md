# <img src="https://github.com/carlosborca/CrystaLattE/blob/master/media/logo/Logo.png" height=260>

[//]: # CrystaLattE
[//]: # ==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/carlosborca/CrystaLattE.png)](https://travis-ci.org/carlosborca/CrystaLattE)
[//]: # [![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/REPLACE_WITH_APPVEYOR_LINK/branch/master?svg=true)](https://ci.appveyor.com/project/REPLACE_WITH_OWNER_ACCOUNT/CrystaLattE/branch/master)
[![codecov](https://codecov.io/gh/carlosborca/CrystaLattE/branch/master/graph/badge.svg)](https://codecov.io/gh/carlosborca/CrystaLattE/branch/master)

Automated calculation of crystal lattice energies with the many-body expansion

## Overview

This is the overview.

## Installation

Minimal set of commands to install and run CrystaLattE on Linux or MacOS:

### Install miniconda:

Get the installer from the Anaconda website:

```wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh```

Check the installer (Optional):

```sha256sum Miniconda3-latest-Linux-x86_64.sh```

Run the installer following the on-screen instructions:

```bash Miniconda3-latest-Linux-x86_64.sh```

After the installation is complete, restart the shell.

Disable automatic activation of the _base_ conda environment _(Optional)_:

```conda config --set auto_activate_base false```

Test the conda installation:

```conda list```

### Create a conda environment for CrystaLattE:

```conda create -n CrystaLattE python=3.7 psi4 psi4-rt pycifrw -c psi4/label/dev -c psi4 -c conda-forge```

Note: The `v2rdm_casscf` plugin of PSI4 can cause errors when trying to import PSI4 from CrystaLattE. In that case, remove it.

```conda remove v2rdm_casscf```

### Activate the CrystaLattE conda environment:

```conda activate CrystaLattE```

### Clone CrystaLattE from its GitHub repository:

```git clone https://github.com/carlosborca/CrystaLattE.git```

### Go to the root directory (`CrystaLattE/`) and run the test suite:

```pytest```



### Copyright

Copyright (c) 2019, Carlos H. Borca


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.0.
