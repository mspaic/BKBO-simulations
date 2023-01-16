Description: 
Slightly updated version of the program "Scatty" written by Joe Paddison
(link: https://joepaddison.com/software/) and published alongside 
the paper http://scripts.iucr.org/cgi-bin/paper?S2053273318015632. 
Scatty is a program for FFT based calculation of diffuse scattering
from atomistic models. 
------------------------------------------------------------------------
Notable modifications: 
1. Performance critical sections of the code were parallelized 
using OpenMP directives.


2. Section of the original code that calculates the products 
of displacements required for Taylor expansion was rewritten 
from scratch using a significantly faster recursive algorithm. 

3. Both terminal and file output formatting was changed where
necessary from generic list directed output to more explicit 
formatting. This was done to reduce the output format dependence 
on the choice of compiler. 

4. Small variable scoping issues preventing the compilation 
using Intel compiler were treated by renaming. 

5. Code was reformatted using fprettify command: 
-> fprettify -i 3 --case 2 2 2 2  scatty.f90
-------------------------------------------------------------------------
-Compile using gfortran: 

gfortran -O3 -fopenmp -march=native singleton.f90 scatty.f90 -o scatty

or simply run "make" in the "programs" directory.
The latter option will explicitly build .o object files
and then link them to build the executable called "scatty". 

-Compile using Intel compiler: 

ifort  -O3 -qopenmp -march=native singleton.f90 scatty.f90 -o scatty

-for more aggressive (but numerically unsafe) optimization -Ofast flag 
can be used with both ifort and gfortran

-If stack overflow occurs for large array allocations 
(manifesting as a runtime segmentation fault) 
add -heap-arrays flag (Intel) or -fno-stack-arrays (gfortran)
to prevent the program from allocating large arrays on stack

--------------------------------------------------------------------------
-For all instructions on running scatty and preparing input files see 
"Instructions.pdf". To be able to run scatty from any point in the 
directory tree OS environment varible needs to be set which 'points'
to the location of the executable. On Linux this can be 
done by appending 

"export PATH=$PATH:/path/to/Installation/directory"

to ~/.bashrc file. 

NOTE: If diffuse scattering is calculated for displacive disorder,
memory overflow can occur if the Taylor expansion is truncated at high 
order to satisfy the maximum error constraint. In that case, maximum
EXPANSION ORDER needs to be specified in addition to EXPANSION_MAX_ERROR
in order to prevent Scatty from allocating too much memory. Scatty will 
automatically compensate by using direct Fourier summation where necessary 
to keep the calculation error below specified value. 
-------------------------------------------------------------------------- 
