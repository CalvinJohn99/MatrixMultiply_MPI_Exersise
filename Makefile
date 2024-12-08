CXX = g++
CXXFLAGS = -g -O3 -Wall -Wextra -Wpedantic

MPICXX?=/usr/bin/mpic++ #/usr/local/bin/mpic++
#MPILIB=-L/usr/local/lib/ -lmpi
#MPIINCLUDE=-I/usr/local/include/

TARGETS= serial mpi

default : all

serial: matrixMultiply.cpp
	$(CXX) $(CXXFLAGS) matrixMultiply.cpp -o $@

mpi: matrixMultiplyMPI.cpp
	$(MPICXX) $(CXXFLAGS) $(MPILIB) $(MPIINCLUDE) matrixMultiplyMPI.cpp -o $@


all : $(TARGETS)

clean:
	rm -f $(TARGETS) *.o *.exe
