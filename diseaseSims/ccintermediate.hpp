#ifndef __CCINTERMEDIATE_HPP__
#define __CCINTERMEDIATE_HPP__

#include <vector>
#include <list>
#include <functional>
#include <iosfwd>
#include <Sequence/SimData.hpp>
#include <diseaseSims/mutation_with_age.hpp>
#include <zlib.h>

struct cc_intermediate
{
  unsigned ncontrols,ncases;
  Sequence::SimData neutral,causative; //genotypes
  std::vector<char> min_n,min_c; //minor alleles
  std::vector<double> G,E; //genetic and random bits of phenotype
  std::vector<unsigned> control_ids,case_ids; //the indexes relating individuals back to general pop
  cc_intermediate(void);

  //IO routines
  std::ostream & buffer( std::ostream & o ) const;
};

std::ostream & operator<<(std::ostream &, const cc_intermediate & );

cc_intermediate process_population( const std::vector< std::pair<glist::iterator,glist::iterator> > & diploids,
				    const std::vector<std::pair<double,double> > & phenotypes,
				    const std::vector<unsigned> & put_controls,
				    const std::vector<unsigned> & put_cases,
				    const unsigned & ncontrols,
				    const unsigned & ncases);

std::pair<double,double> phenosums(const std::vector<std::pair<double,double> > & phenos, const double & case_proportion, double * cutoff);

void grab_putative_CC( const std::pair<double,double> & mean_sd,
		      const std::vector<std::pair<double,double> > & phenotypes,
		      const double & crange,
		      const double & cutoff,
		      std::vector<unsigned> & put_controls,
 		      std::vector<unsigned> & put_cases );
#endif
