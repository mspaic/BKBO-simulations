#!/usr/bin/env python
from setuptools import setup, Extension, distutils
import pybind11

compiler_flags = ['-O3','-march=native'] 
          

bvs_module = Extension(
    'BvsOptimizer',
    sources=['BvsOptimizer.cpp'],
    include_dirs=[pybind11.get_include()],
    language='c++',
    extra_compile_args=compiler_flags 
    )

setup(
    name='BvsOptimizer',
    author='Marin Spaic',
    version='1.0',
    description="A simple Pybind11 wrapped C++ module for atomistic Monte Carlo modeling"
                   "of local structure of disordered materials. It implements a "
                   "few simple interatomic interaction potentials and "
                   "provides the possibility of Monte Carlo "
                   "optimization of atomic oxidation states using bond valence sum" 
                   "method from crystal chemistry, thereby providing a simple way of" 
                   "simulating the local bonding geometry of atoms within lattices. ", 
    ext_modules=[bvs_module]

)
