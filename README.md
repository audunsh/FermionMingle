## Fermion Mingle - A quantum many body solver
### Audun Skau Hansen | Comp-Phys | UiO | 2014
This code is made for educational purposes as part of my masters thesis in computational physics. While most of the functionality is thoroughly tested and benchmarked against similar software, there are some major performance limitations when calculating large systems or basis sets. Some optimalization is planned to be implemented in the nearest future.

I have named the project "Fermion Mingle", as the code enables calculation of fermionic interactions.

Parts of the code was written in collaboration with Gøran Brekke Svaland.

Currently the code has the following functionality:
- Gaussian basis sets and integral calculations
- Restricted Hartree-Fock
- Unrestricted Hartree-Fock
- Coupled Cluster (CCD, CCSD)
- Electron density evaluation

Some functions are implemented in python:
- CCalgebra - An environment for coupled cluster calculations and code generation
- TurbomoleConverter - A script for generating C++ code from turbomole files

