#ifndef __EYRE_WALKER_FITNESS_MODEL_HPP__
#define __EYRE_WALKER_FITNESS_MODEL_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <utility>
#include <algorithm>
#include <numeric>

struct EWfitness 
{
  typedef double result_type;
  template< typename gamete_type,typename mcont_t >
  inline double operator()(const gamete_type & g1, const gamete_type & g2,
			   const mcont_t & mutations) const
  {
    double sum = std::accumulate(g1.smutations.begin(),
				 g1.smutations.end(),
				 0.,
				 [&mutations]( const double & a,
					       const std::size_t i ) { return a + mutations[i].s; } )
      + std::accumulate(g2.smutations.begin(),
			g2.smutations.end(),
			0.,
			[&mutations]( const double & a,
				      const std::size_t i ) { return a + mutations[i].s; } );
    return std::max(0.,1.-sum);
  }
};

#endif
