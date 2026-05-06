# Makefile
# Author: KK

CXXFLAGS = -std=c++20 -g -O2\
	-fopenmp -Wpedantic -Wall -Wextra -Werror\
	-Wdisabled-optimization\
	-Wlogical-op\
	-Wmissing-declarations\
	-Wmissing-include-dirs\
	-Wredundant-decls\
	-Wshadow\
	-Wswitch-default\
	-Wsign-conversion\
	-Wfloat-conversion
CXX = g++ # icpc
OBJS = src/WignerSolver.o \
	src/WignerIO.o \
	src/WignerTools.o \
	src/Poisson1D.o \
	src/TestPoisson1D.o \
	main.o
LDLIBS = -larmadillo -lsuperlu -lopenblas -lm -fopenmp
#-lblas <-> -lopenblas
#-lsuperlu -larmadillo -lopenblas -lm -fopenmp #-m64 -I${MKLROOT}/include # -mkl
#-lsuperlu -L/home/karol/intel/compilers_and_libraries_2019.0.117/linux/mkl/lib -llapack -L/opt/OpenBLAS/lib/ -lm  -lopenblas -lm
# -lm -- math library

run.out: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(OBJS) run.out
