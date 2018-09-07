## using openmpi

### try to use Intel MPI

### use openmpi
# path_openmpi=/opt/openmpi-1.8.1/bin
# CXXFLAGS=-Wall -fopenmp -O3

### use intelmpi
CXXFLAGS=-Wall -qopenmp -O3 

# CXXFLAGSO1=-Wall -std=gnu++0x -fopenmp -O1
# CXXFLAGSO3=-Wall -std=gnu++0x -fopenmp -O3
# CXXDEBUGFLAGS=-Wall -std=gnu++0x -fopenmp -g

all:
	# ${path_openmpi}/mpicxx ${CXXFLAGS} -o fascia fascia-mpi.cpp
	mpiicc ${CXXFLAGS} -o fascia fascia-mpi.cpp
	# mpicxx ${CXXFLAGS} -o fascia fascia-mpi.cpp

clean:
	rm -f fascia 
