#ifndef __GENE_BASED_MODEL_HPP__
#define __GENE_BASED_MODEL_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <utility>

//Functions for the Thornton, Foran, and Long (2013) model
struct disease_effect
{
  typedef double result_type;
  template< typename iterator_type >
  inline std::pair<double,double> operator()(const iterator_type & g1, const iterator_type & g2,
					     const double & sd, gsl_rng * r) const
  {
    //The effect of each allele is additive across mutations
    double e1 = 0.,e2=0.;
    typename  iterator_type::value_type::mutation_container::const_iterator itr;
    for(itr = g1->smutations.begin() ; itr != g1->smutations.end() ; ++itr)
      {
	e1 += (*itr)->s;
      }
    for(itr = g2->smutations.begin() ; itr != g2->smutations.end() ; ++itr)
      {
	e2 += (*itr)->s;
      }
    double effect = pow( e1*e2, 0.5 );
    return std::make_pair(effect, gsl_ran_gaussian(r,sd));
  }
};

//calculates the fitess of a diploid
struct disease_effect_to_fitness
{
  typedef double result_type;
  template< typename iterator_type >
  inline double operator()(const iterator_type & g1, const iterator_type & g2,
			   const double & sd, const double & sd_s,gsl_rng * r) const
  {
    std::pair<double,double> effect = disease_effect()(g1,g2,sd,r);
    double fitness = exp( (-1. * pow(effect.first+effect.second,2.))/(2.*pow(sd_s,2)) );
    return ( fitness );
  }
};

#endif
