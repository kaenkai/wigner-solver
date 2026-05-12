# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++20 -g -O2 -Isrc -fopenmp\
	-Wpedantic -Wall -Wextra -Werror\
	-Wdisabled-optimization\
	-Wlogical-op\
	-Wmissing-declarations\
	-Wmissing-include-dirs\
	-Wredundant-decls\
	-Wshadow\
	-Wswitch-default\
	-Wsign-conversion\
	-Wfloat-conversion

# Object files
OBJS = src/WignerSolver.o src/WignerIO.o src/WignerTools.o src/Poisson1D.o\
	src/test/TestDerivative.o src/test/TestPoisson1D.o src/test/TestWignerSolver.o main.o

# Libraries
LIBS = -larmadillo -lsuperlu -lopenblas -lm

run: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f $(OBJS) run
