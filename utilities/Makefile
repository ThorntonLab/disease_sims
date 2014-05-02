CXX=c++
DEBUG=
CXXFLAGS=-I. -I.. -O2 -Wall -W $(DEBUG)

all: phenoburden.o
	$(CXX) $(CXXFLAGS) -o phenoburden phenoburden.o -lboost_system -lboost_program_options $(LDFLAGS)

phenoburden.o: ../readCC.hpp
