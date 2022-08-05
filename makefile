SOURCES=non_blocking_n_body_barnes_hut.cpp
# the compiler: mpic++ for mpi
MPICXX = mpic++

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
MPICXX_FLAGS  = -std=c++14 -march=native -g -Wall -pedantic -Wall -Wfatal-errors -Wextra -Wno-unused-parameter -Wno-unused-variable
OPTIM_FLAGS = -fopenmp -fopenmp-simd -Ofast 
HDFQL_FLAGS = -I./hdfql-2.2.0/include ./hdfql-2.2.0/wrapper/cpp/libHDFql.a  -lmpi -ldl
HDFQL_FLAGS_NOMPI = -I./hdfql-2.2.0/include ./hdfql-2.2.0/wrapper/cpp/libHDFql.a  -lmpi -ldl

.PHONY: default mpionly clean
default:
	$(MPICXX) $(SOURCES) $(OPTIM_FLAGS) $(MPICXX_FLAGS) $(HDFQL_FLAGS) -o nbody_all.o

# Doesnt work
mpionly:
	$(MPICXX) $(SOURCES) $(MPICXX_FLAGS) -o nbody_mpi.o

clean:
	rm *.o
