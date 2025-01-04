# WignerSolver ver. `0.250104`

Program for solving Boltzmann-Poisson and Wigner-Poisson set of equations.

* Written in C++.
* Uses armadillo C++ library for linear algebra & scientific computing. 
* Due to the use of spare matrices SuperLU library is needed.
* Classes necessary for computation are located in 'src' folder.
* Program is compiled through makefile.
* To disable multi-threading, delete '-fopenmp' in makefile.
* Google drive backup: 28/12/2024.

## Compilation

Compile the program using the following command: `make`

To delete the object files and the executable, use the following command: `make clean`

## Usage

To run the program, use the following command: `./run.out`

The inital parameters are set in the `main.cpp` file.

## Dependencies

The program uses the following libraries:
* Armadillo C++ library for linear algebra & scientific computing.
* SuperLU library for solving sparse linear systems.