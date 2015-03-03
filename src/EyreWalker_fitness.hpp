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
  template< typename iterator_type >
  inline double operator()(const iterator_type & g1, const iterator_type & g2) const
  {
    using itr = typename  iterator_type::value_type::mutation_container::const_iterator::value_type;
    double sum = std::accumulate(g1->smutations.begin(),
				 g1->smutations.end(),
				 0.,
				 []( const double & a,
				     const itr & i ) { return a + i->s; } )
      + std::accumulate(g2->smutations.begin(),
			g2->smutations.end(),
			0.,
			[]( const double & a,
			    const itr & i ) { return a + i->s; } );
    return std::max(0.,sum);
  }
};

#endif
