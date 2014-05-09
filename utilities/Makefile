CXX=c++
DEBUG=
CXXFLAGS=-I. -I.. -I../common -O2 -Wall -W $(DEBUG)

all: phenoburden.o CCburden.o ../readCC.o sfs.o sfsSample.o
	$(CXX) $(CXXFLAGS) -o phenoburden phenoburden.o ../locking_routines.o $(LDFLAGS) -lboost_system -lboost_program_options 
	$(CXX) $(CXXFLAGS) -o CCburden CCburden.o ../readCC.o ../locking_routines.o $(LDFLAGS) -lboost_system -lboost_program_options 
	$(CXX) $(CXXFLAGS) -o sfs sfs.o ../common/readSimOutput.o $(LDFLAGS) -lboost_system -lboost_program_options 
	$(CXX) $(CXXFLAGS) -o sfsSample sfsSample.o $(LDFLAGS) -lsequence -lz -lboost_system -lboost_program_options -lgsl -lgslcblas

phenoburden.o: ../locking_routines.hpp ../mutation_with_age.hpp
CCburden.o: ../readCC.o ../locking_routines.hpp ../mutation_with_age.hpp
sfs.o: ../common/readSimOutput.hpp ../common/readSimOutput.o
sfsSample.o: ../mutation_with_age.hpp