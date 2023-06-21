
## Description
- Compile Ising.cpp (e.g. with GNU compiler):

   `g++ -Wall -O3 -march=native -fopenmp -o Ising.out Ising.cpp`

- Run with `./Ising.out` and enter simulation parameters when prompted. Option `get_samples`
  generates a specified number of equilibrium configurations while the option `sweep_doping`
  is used to measure the doping dependence of the order parameter (staggered magnetization).

  Note: Supercell size is fixed by setting the global variable "size" in Ising.cpp.
 
- breathing_distortion.py imports configurations created by Ising.out
  and creates an accompanying octahedral "breathing" distortion using a
  Bi-O harmonic interaction model.

  Note: Javelin does not distinguish oxidation states. Therefore, we represent
  Bi5+ and Bi3+ by atoms with similar form factors.

  Bi5+ -> s = +1  (composition y = (1+x)/2) -> Represented in javelin by Pt

  Bi3+ -> s = -1  (composition 1-y) -> Represented in javelin by Hg

  K -> sigma = 1  (composition x)

  Ba-> sigma = 0  (composition 1-x)

  Couplings:

  J1 > 0  (Bi3+ and Bi5+ anticorrelated -> CDW formation)

  J2 < 0  (K locally favours the Bi5+ state)

- ising.ipynb contains functions used for analysis of short-range order
  created by ising model simulation

- To calculate diffuse scattering in the L=1/2 plane use diffuse.ipynb to import
  structures with relaxed Bi-O bonds and generate input files for scatty (FFT based
  diffuse scattering calculator) -> see Dependencies/Scatty subfolder for instructions
  on using scatty




