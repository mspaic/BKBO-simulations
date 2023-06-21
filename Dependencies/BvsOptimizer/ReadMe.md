# BvsOptimizer

## Description

A simple Pybind11 wrapped C++ module for atomistic Monte Carlo
modeling  of local structure of disordered materials. It
implements a few simple interatomic interaction potentials
and provides the possibility of Monte optimization of atomic
oxidation states using bond valence sum method from crystal
chemistry, thereby providing a simple way of simulating the
local bonding geometry of atoms within lattices.

## Build and install
 - Install dependencies listed in environment.yml file using e.g. conda

   `conda env update --file environment.yml`

 - Build with setuptools

   `python setup.py build_ext --inplace`

 - Install with pip

   `pip install -e .`


NOTE: Compilation flags in setup.py presuppose that GNU compilers
are set as default C/C++ compilers. If other compilers (msvc,
clang, icpc) are used, compiler flags will need to be modified
accordingly.



