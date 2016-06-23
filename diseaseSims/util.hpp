#ifndef __DISEASESIMS_UTIL_HPP__
#define __DISEASESIMS_UTIL_HPP__

#include <functional>
#include <utility>
#include <vector>
#include <cstdlib>
#include <zlib.h>
#include <diseaseSims/mutation_with_age.hpp>

using Gfxn_t = std::function<double(const std::size_t first,
				    const std::size_t second,
				    const poptype::gcont_t & gametes,
				    const poptype::mcont_t & mutations)>;
//pair = index of mutation in pop.mutations, count of that mutation in a diploid
using vmcount_t = std::vector<std::pair<std::size_t,std::int8_t> >;

poptype readPop( gzFile gzin );
std::vector<double> getG( const poptype & pop,
			  const Gfxn_t & dipG );
//Counts of mutations per diploid
vmcount_t get_mut_counts( const poptype & pop,
			  const std::size_t dipindex,
			  const bool & selected = true);

#endif
