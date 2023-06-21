# Scatty
## Description
Modified version of the program [Scatty](https://joepaddison.com/software/) written by Joe Paddison and published alongside the paper
[Paddison, J. A. M. (2019)](http://scripts.iucr.org/cgi-bin/paper?S2053273318015632).
Scatty is a program for FFT based calculation of diffuse scattering from atomistic models.

## Notable modifications
- Performance critical sections of the code were parallelized
  using OpenMP directives.

- Section of the original code that calculates the products
  of displacements required for Taylor expansion was rewritten
  from scratch using a significantly faster recursive algorithm.

- Both terminal and file output formatting was changed where
  necessary from generic list directed output to more explicit
  formatting. This was done to reduce the output format dependence
  on the choice of compiler and to make the ouput more readable.

- Small variable scoping issues preventing the compilation
  using Intel compiler were treated by renaming.

- Code was reformatted using fprettify command:

  `fprettify -i 3 --case 2 2 2 2  scatty.f90`

## Bulid
 - Using gfortran:

    `gfortran -O3 -fopenmp -march=native singleton.f90 scatty.f90 -o scatty`

    or simply run `make` in the "programs" directory.
    The latter option will explicitly build .o object files
    and then link them to build the executable named "scatty". 

- Using Intel compiler:

    `ifort  -O3 -qopenmp -march=native singleton.f90 scatty.f90 -o scatty`


## Usage
   - For all instructions on running scatty and preparing input files see
    "Instructions.pdf". To be able to run scatty from any point in the 
    directory tree, OS environment variable needs to be set which 'points'
    to the location of the executable. On Linux this can be 
    done by appending 

      `export PATH=$PATH:/path/to/Installation/directory`

      to `~/.bashrc` file.

## Known issues

  - Segmentation fault sometimes occurs at the very end of execution for no
    obvious reason when scatty is compiled with Intel compiler. Usually this occurs
    when large arrays are passed to subroutine 'write fit' and is possibly
    due to creation of large temporary arrays on stack during function call.
    A simple solution is to add -heap-arrays flag (Intel) to prevent the program
    from using stack for temporary arrays. This solves the problem, but makes
    the compilation painfully slow. On the other hand, performance penalty
    for doing this is not significant.


