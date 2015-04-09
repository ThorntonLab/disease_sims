#ifndef __DISEASESIMS_UTIL_HPP__
#define __DISEASESIMS_UTIL_HPP__

#include <functional>
#include <utility>
#include <vector>
#include <cstdlib>
#include <diseaseSims/mutation_with_age.hpp>

using Gfxn_t = std::function<double(const glist::const_iterator &,
				    const glist::const_iterator &)>;
using vmcount_t = std::vector<std::pair<mlist::iterator,std::int8_t> >;

struct popstruct
{
  mlist mutations;
  glist gametes;
  dipvector diploids;
};

popstruct readPop( gzFile gzin );
std::vector<double> getG( const dipvector & diploids,
		     const Gfxn_t & dipG );
//Counts of mutations per diploid
vmcount_t get_mut_counts( const glist::const_iterator & g1,
			  const glist::const_iterator & g2,
			  const bool & selected = true);

#endif
