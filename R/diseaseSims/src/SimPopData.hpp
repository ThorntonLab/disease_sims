#ifndef __SIMPOPDATA_HPP__
#define __SIMPOPDATA_HPP__

#include <diseaseSims/mutation_with_age.hpp>

class SimPopData
{
private:
  void readPopWrapper(const char * filename,
		      const unsigned long & offset);
public:
  poptype pop;

  SimPopData(const std::string & filename,
	     const unsigned long & offset);

  Rcpp::List neutralGenotype(const size_t & i) const;
  Rcpp::List selectedGenotype(const size_t & i) const;
  poptype::dipvector_t::size_type popsize() const;
};

Rcpp::NumericMatrix getPheno(const char * filename,
			     const unsigned long & offset);

#endif
