#ifndef __CCINTERMEDIATE_HPP__
#define __CCINTERMEDIATE_HPP__

#include <vector>
#include <list>
#include <functional>
#include <iosfwd>
#include <Sequence/SimData.hpp>
#include <mutation_with_age.hpp>
#include <zlib.h>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

struct cc_intermediate
{
  //allow intrusive serialization
  friend class boost::serialization::access;
  unsigned ncontrols,ncases;
  Sequence::SimData neutral,causative; //genotypes
  std::vector<char> min_n,min_c; //minor alleles
  std::vector<double> G,E; //genetic and random bits of phenotype
  std::vector<unsigned> control_ids,case_ids; //the indexes relating individuals back to general pop
  cc_intermediate(void);

  //IO routines
  std::ostream & buffer( std::ostream & o ) const;
  
  //boost::serialization -- will help us converting data in R, I hope.
  template<typename Archive>
  void serialize( Archive & ar, const unsigned int version )
  {
    ar & ncontrols;
    ar & ncases;
    //As of libseq 1.8.4, PolyTable inherits pair< vector<double>, vector<string>, allowing trivial serialization
    ar & neutral.first;
    ar & neutral.second;
    ar & causative.first;
    ar & causative.second;
    ar & min_n;
    ar & min_c;
    ar & G;
    ar & E;
    ar & control_ids;
    ar & case_ids;
  }
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
