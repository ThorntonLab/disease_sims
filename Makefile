CXX=c++
CXXFLAGS=-O2 -Wall -W -I.

TARGETS=TFL2013_ind

all: TFL2013.o TFL2013_ind.o
	$(CXX) $(CXXFLAGS) -o TFL2013_ind TFL2013_ind.o -lboost_system -lboost_program_options -lgsl -lgslcblas $(LDFLAGS)
clean: 
	rm -f *.o $(TARGET)

TFL2013_ind.o: mutation_with_age.hpp TFL_fitness_models.hpp gene_based_model.hpp multiplicative_model.hpp