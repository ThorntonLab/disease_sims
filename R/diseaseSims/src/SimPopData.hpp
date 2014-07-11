#ifndef __SIMPOPDATA_HPP__
#define __SIMPOPDATA_HPP__

#include <mutation_with_age.hpp>

class SimPopData
{
private:
  void readPop(const char * filename,
	       const unsigned long & offset);
public:
  mlist mutations;
  glist gametes;
  typedef boost::container::vector< std::pair<glist::iterator,glist::iterator> > dips;
  dips diploids;

  SimPopData(const std::string & filename,
	     const unsigned long & offset) : mutations(mlist()),
					     gametes(glist()),
					     diploids(dips())
  {
    this->readPop(filename.c_str(),offset);
  }

  Rcpp::List neutralGenotype(const size_t & i) const;
  Rcpp::List selectedGenotype(const size_t & i) const;
  dips::size_type popsize() const;
};

#endif
