#ifndef __CCINTERMEDIATE_HPP__
#define __CCINTERMEDIATE_HPP__

#include <vector>
#include <list>
#include <functional>

#include <Sequence/SimData.hpp>
#include <mutation_with_age.hpp>

struct cc_intermediate
{
  Sequence::SimData neutral,causative;
  std::vector<char> min_n,min_c; //minor alleles
  std::vector<std::pair<double,double> > phenotypes;
  
  cc_intermediate(void);
};

cc_intermediate process_population( const std::vector< std::pair<glist::iterator,glist::iterator> > & diploids,
				    const std::vector<std::pair<double,double> > & phenotypes,
				    const std::vector<size_t> & put_controls,
				    const std::vector<size_t> & put_cases,
				    const unsigned & ncontrols,
				    const unsigned & ncases);

#endif
