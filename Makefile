CXX=c++
CXXFLAGS=-O2 -Wall -W -I.

TARGETS=TFL2013_ind

all: TFL2013.o TFL2013_ind.o make_case_control.o readCC.o esm_chisq_zscore.o esm.o chisq_per_marker.o locking_routines.o esm_filter_sites.o 
	$(CXX) $(CXXFLAGS) -o TFL2013_ind TFL2013_ind.o -lboost_system -lboost_program_options -lgsl -lgslcblas $(LDFLAGS)
clean: 
	rm -f *.o $(TARGETS)

TFL2013_ind.o: mutation_with_age.hpp TFL_fitness_models.hpp gene_based_model.hpp multiplicative_model.hpp

make_case_control.o: mutation_with_age.hpp

readCC.o: readCC.hpp

esm_chisq_zscore.o: readCC.hpp

chisq_per_marker.o: chisq_per_marker.hpp

esm.o: esm.hpp

esm_filter_sites.o: esm_filter_sites.hpp

locking_routines.o: locking_routines.hpp