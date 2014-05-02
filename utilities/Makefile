CXX=c++
DEBUG=
CXXFLAGS=-I. -I.. -O2 -Wall -W $(DEBUG)

all: phenoburden.o CCburden.o ../readCC.o
	$(CXX) $(CXXFLAGS) -o phenoburden phenoburden.o ../locking_routines.o -lboost_system -lboost_program_options $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o CCburden CCburden.o ../readCC.o ../locking_routines.o -lboost_system -lboost_program_options $(LDFLAGS)

phenoburden.o: ../locking_routines.hpp ../mutation_with_age.hpp
CCburden.o: ../readCC.o ../locking_routines.hpp
