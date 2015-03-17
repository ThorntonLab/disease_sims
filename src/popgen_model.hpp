#ifndef __POPGEN_MODEL_HPP__
#define __POPGEN_MODEL_HPP__

#include <diseaseSims/traitValues.hpp>

struct popgen_disease_effect_to_fitness
{
  typedef double result_type;
  template< typename iterator_type >
  inline double operator()(const iterator_type & g1, const iterator_type & g2,
			   const double & dominance,
			   const double & sd, const double & sd_s,const double & optimum,gsl_rng * r) const
  {
    double genetic = popgen_phenotype()(g1,g2,dominance);
    double noise = (sd>0.)?gsl_ran_gaussian(r,sd):0.;
    //Subtract 1 so that phenotype has mean 1 and std_dev sd_s
    double fitness = exp( (-1. * pow(std::abs(genetic+noise)-1.-optimum,2.))/(2.*pow(sd_s,2)) );
    return ( fitness );
  }
};

#endif
