CXX=c++
OPT=-O2
DEBUG=-DNDEBUG
#DEBUG=
CXXFLAGS=$(OPT) $(DEBUG) -Wall -W -I.

TARGETS=TFL2013_ind make_case_control esm_chisq_zscore atomic_locker/atomic_locker

all: TFL2013.o TFL2013_ind.o make_case_control.o readCC.o esm_chisq_zscore.o esm.o chisq_per_marker.o locking_routines.o esm_filter_sites.o ccintermediate.o
	$(CXX) $(CXXFLAGS) -o TFL2013 TFL2013.o -lboost_system -lboost_iostreams -lgsl -lgslcblas -lsequence $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o TFL2013_ind TFL2013_ind.o -lboost_system -lboost_program_options -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o make_case_control make_case_control.o ccintermediate.o locking_routines.o -lsequence -lz -lboost_system -lboost_program_options -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CXXFLAGS) -o esm_chisq_zscore esm_chisq_zscore.o esm.o chisq_per_marker.o esm_filter_sites.o readCC.o locking_routines.o -lboost_program_options -lboost_system -lboost_thread -lgsl -lgslcblas $(LDFLAGS)
	git submodule init
	git submodule update
	cd atomic_locker && make
clean: 
	rm -f *.o atomic_locker/*.o $(TARGETS)

TFL2013_ind.o: mutation_with_age.hpp TFL_fitness_models.hpp gene_based_model.hpp multiplicative_model.hpp

make_case_control.o: mutation_with_age.hpp ccintermediate.hpp

readCC.o: readCC.hpp

esm_chisq_zscore.o: readCC.hpp gwas_stats.hpp locking_routines.hpp

chisq_per_marker.o: chisq_per_marker.hpp

esm.o: esm.hpp

esm_filter_sites.o: esm_filter_sites.hpp

locking_routines.o: locking_routines.hpp