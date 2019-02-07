# calc2DIR

## What is it?
This program calculates linear and 2D spectra for *coupled* chromophores.
The calculation is based on material found in "Concepts and methods of 2D Infrared Spectroscopy" by Peter Hamm and Martin Zanni and "Principles of Nonlinear Optical Spectroscopy" by Mukamel.

## How do I use it
The program can be built using the supplied Makefile.
The [MKL](https://software.intel.com/en-us/mkl/choose-download/linux) and [fftw3](http://www.fftw.org/download.html) libraries are required.

The program can be run from the command line using the command:

calc2DIR.exe input.inp

where input.inp is an input file.
An example input file is included in this repository.
The program requires energy and dipole binary files, the names of which are specified in the input file.
These files should be created by the user.
For each frame in the trajectory, the energy file should contain first the number of the current frame as an integer, followed by the one-exciton hamiltonian in upper tridiagonal format.
The two-exciton Hamiltonian is built automatically in the program, given the anharmonicity supplied in the input file.
For each frame in the trajectory, the dipole file should contain first the number of the current frame as an integer, followed by the x, y and z components of the chromophore's ground-to-one exciton dipole moment as floats.
The one-to-two exciton dipole moments are generated automatically in the program assuming they behave like in a harmonic oscillator.
If it is not clear how these files should be formated from this README, check the source code subroutines readDfile and readEfile to see how the program handles reading these files.

Also included is a tutorial folder.
The run script will generate energy and dipole files according to the parameters specified in that file.
The back end program stochastic generates correlated random energies and writes energy and dipole binary files that can be used as input.
This program was written by T.L.C Jansen to simulate the Hamiltonian for a coupled dimer.
The repository for his 2D IR code can be found [here](https://github.com/GHlacour/NISE_2015).

## Contact
This README is very brief. If questions arise, please contact me at nicholasjhestand at gmail.com.
