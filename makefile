# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++20 -g -O2 -Isrc -fopenmp -Werror -Wpedantic -Wall -Wextra

# Object files
OBJS = src/WignerSolver.o src/WignerIO.o src/WignerTools.o src/Poisson1D.o src/test/TestPoisson1D.o src/test/TestWignerSolver.o main.o

# Libraries
LIBS = -larmadillo -lsuperlu -lopenblas -lm

run: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f $(OBJS) run
