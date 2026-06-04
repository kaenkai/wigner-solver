# wigner-solver

Program for solving Boltzmann-Poisson (classical) and Wigner-Poisson (quantum) electron transport systems.

* Written in C++, and employs **Armadillo C++ library for linear algebra & scientific computing**. 
* Due to the use of spare matrices **SuperLU** library is needed.
* Classes necessary for computation are located in *src* folder.
* Program is compiled through makefile.
* To disable multi-threading, delete `-fopenmp` in makefile.

## Compilation

* Compile the program using the following command: `make`.
* To delete the object files and the executable, use the following command: `make clean`.

## Usage

* To run the program, use the following command: `./run`.
* The system parameters are set throug the `main.cpp` file.

## Dependencies

* Armadillo C++ library for linear algebra & scientific computing.
* SuperLU library for solving sparse linear systems.

## TODO

* Fix stability issues in Boltzmann-Poisson system
* Clean up unused code
