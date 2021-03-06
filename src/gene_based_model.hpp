#ifndef __GENE_BASED_MODEL_HPP__
#define __GENE_BASED_MODEL_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <utility>
#include <algorithm>

//Functions for the Thornton, Foran, and Long (2013) model
struct TFL2013_recessive_disease_effect
{
  typedef double result_type;
  template< typename iterator_type >
  inline std::pair<double,double> operator()(const iterator_type & g1, const iterator_type & g2,
					     const double & sd, gsl_rng * r) const
  {
    //The effect of each allele is additive across mutations
    double e1 = 0.,e2=0.;
    using itr = typename  iterator_type::value_type::mutation_container::const_iterator::value_type;
    std::for_each(g1->smutations.cbegin(),g1->smutations.cend(),
		  [&](const itr & __i) { e1 += __i->s; });
    std::for_each(g2->smutations.cbegin(),g2->smutations.cend(),
		  [&](const itr & __i) { e2 += __i->s; });
    double effect = pow( e1*e2, 0.5 );
    return std::make_pair(effect, (sd>0.)?gsl_ran_gaussian(r,sd):0.);
  }
};

struct TFL2013_additive_disease_effect
{
  typedef double result_type;
  template< typename iterator_type >
  inline std::pair<double,double> operator()(const iterator_type & g1, const iterator_type & g2,
					     const double & sd, gsl_rng * r) const
  {
    //The effect of each allele is additive across mutations
    double e1 = 0.,e2=0.;
    using itr = typename  iterator_type::value_type::mutation_container::const_iterator::value_type;
    std::for_each(g1->smutations.cbegin(),g1->smutations.cend(),
		  [&](const itr & __i) { e1 += __i->s; });
    std::for_each(g2->smutations.cbegin(),g2->smutations.cend(),
		  [&](const itr & __i) { e2 += __i->s; });
    double effect = e1 + e2;
    return std::make_pair(effect, (sd>0.)?gsl_ran_gaussian(r,sd):0.);
  }
};

//calculates the fitess of a diploid under the recessive model of TFL (2013)
struct TFL2013_recessive
{
  typedef double result_type;
  template< typename iterator_type >
  inline double operator()(const iterator_type & g1, const iterator_type & g2,
			   const double & sd, const double & sd_s,const double & optimum,gsl_rng * r) const
  {
    std::pair<double,double> effect = TFL2013_recessive_disease_effect()(g1,g2,sd,r);
    double fitness = exp( (-1. * pow(effect.first+effect.second-optimum,2.))/(2.*pow(sd_s,2)) );
    return ( fitness );
  }
};

/*
  Haplotype effect sizes additive across causative mutations.
  Genetic value of a diploid is additive for the two haplotypes.
 */
struct TFL2013_additive
{
  typedef double result_type;
  template<typename iterator_type>
  inline double operator()(const iterator_type & g1, const iterator_type & g2,
			   const double & sd, const double & sd_s,const double & optimum,gsl_rng * r) const
  {
    std::pair<double,double> effect = TFL2013_additive_disease_effect()(g1,g2,sd,r);
    double fitness = exp( (-1. * pow(effect.first+effect.second-optimum,2.))/(2.*pow(sd_s,2)) );
    return ( fitness );
  }
};



#endif
