[![Build Status](https://travis-ci.com/temken/XXX.svg?token=CWyAeZfiHMD8t4eitDid&branch=master)](https://travis-ci.com/temken/XXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)


# Dirac vs Majorana dark matter

<!-- [![DOI](https://zenodo.org/badge/XXXXXXXX.svg)](https://zenodo.org/badge/latestdoi/XXXXXXXX) -->
<!-- [![arXiv](https://img.shields.io/badge/arXiv-2003.xxxxx-B31B1B.svg)](https://arxiv.org/abs/2003.xxxxx) -->

<!-- <img src="url" width="800"> -->

## GENERAL NOTES

Direct detection experiments which look for sub-GeV dark matter (DM) particles often search for DM-induced electronic transitions inside a target. Assuming one of these experiments would succeed, the next question would be to study the properties of DM.

One question we could ask is if DM particles are their own anti-particles, a so-called Majorana fermion. This code determines the statistical significance with which a successful electron scattering experiment could reject the Majorana hypothesis (using the likelihood ratio test) in favour of the hypothesis of Dirac DM. We assume that the DM interacts with the photon via higher-order electromagnetic moments.

> :warning: **Warning**: In order for this code to produce results, the */data/* folder needs to contain the tabulated atomic response functions, which can be computed with the [DarkARC tool](https://github.com/temken/DarkARC). Without these tables, it is not possible to compute ionization spectra and predictions for signal event rates.

## CONTENT

- */bin/*: Contains the executable.
- */build/*: Will contain the object files after compilation.
- */data/*: The folder for the tabulated atomic response functions.
- */include/*: All the header files can be found here.
- */results/*: Resulting tables and files are saved to this folder.
- */src/*: Contains the source code.

## Installation:

The code can be compiled using the makefile. It might be necessary to adjust the compiler lines and the path to the libraries:

```
#Compiler and compiler flags
CXX := g++
CXXFLAGS := -Wall -std=c++11 
LIB := 
INC := -I include
(...)
```

The code is compiled by running 
```
make
```
from the root directory in the terminal to compile DiracVsMajorana.

Running
```
make clean
```
deletes all object files and executables.


## CITING THIS CODE

If you decide to use this code, please cite the latest archived version,

> [[DOI:10.5281/zenodo.XXXXXXXX]](https://doi.org/10.5281/zenodo.XXXXXXXX)

as well as the original publications,

>Catena, R., Emken, T. , Ravanis J., **Rejecting the Majorana nature of dark matter with electron scattering experiments**, [[arXiv:2003.xxxxx]](https://arxiv.org/abs/2003.xxxxx).

## VERSIONS

- **v1.0** (06/03/2020): Version released with v1 of the preprint [[arXiv:2003.xxxxxv1]](https://arxiv.org/abs/2003.xxxxxv1).

## AUTHORS & CONTACT

The author of this code is Timon Emken.

For questions, bug reports or other suggestions please contact [emken@chalmers.se](mailto:emken@chalmers.se).


## LICENCE

This project is licensed under the MIT License - see the LICENSE file.