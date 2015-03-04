#ifndef __POPGEN_MODEL_HPP__
#define __POPGEN_MODEL_HPP__

#include <fwdpp/fitness_models.hpp>

//Popgen-like model for trait values

struct popgen_phenotype_updater_hom
{
  typedef void result_type;
  template<typename iterator_type>
  inline void operator()(double & fitness, const iterator_type & m1) const
  {
    /*
      Bug fix March 4, 2016
      Old code:     fitness *= (1. + m1->s);
      It was a bug b/c the trait values
      didn't scale as they would for homozygote
      in Risch model.
    */
    fitness *= (1. + 2.*m1->s);
  }
};

struct popgen_phenotype_updater_het
{
  typedef void result_type;
  template<typename iterator_type>
  inline void operator()(double & fitness, const iterator_type & m1,const double & dominance) const
  {
    fitness *= ( 1. + dominance*m1->s );
  }
};

struct popgen_phenotype
{
  typedef double result_type;
  template< typename iterator_type>
  inline double operator()(const iterator_type & g1, const iterator_type & g2,
			   const double & dominance) const
  {
    return KTfwd::site_dependent_fitness()(g1,g2,
					   std::bind(popgen_phenotype_updater_hom(),std::placeholders::_1,std::placeholders::_2),
					   std::bind(popgen_phenotype_updater_het(),std::placeholders::_1,std::placeholders::_2,dominance),
					   1.);
  }
};

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
