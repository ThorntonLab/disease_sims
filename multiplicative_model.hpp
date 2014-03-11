#ifndef __MULTIPLICATIVE_MODEL_HPP__
#define __MULTIPLICATIVE_MODEL_HPP__

#include <fwdpp/fitness_models.hpp>

//Define a model where phenotype is 1, 1+s, (1+s)^2 for a site, then total value is product over sites
struct multiplicative_phenotype_updater_hom
{
  typedef void result_type;
  template<typename iterator_type>
  inline void operator()(double & fitness, const iterator_type & m1) const
  {
    fitness *= ( std::pow(1. + m1->s,2.) );
  }
};

struct multiplicative_phenotype_updater_het
{
  typedef void result_type;
  template<typename iterator_type>
  inline void operator()(double & fitness, const iterator_type & m1) const
  {
    fitness *= ( 1. + m1->s );
  }
};

struct multiplicative_phenotype
{
  typedef double result_type;
  template< typename iterator_type>
  inline double operator()(const iterator_type & g1, const iterator_type & g2) const
  {
    return KTfwd::site_dependent_fitness()(g1,g2,
					   boost::bind(multiplicative_phenotype_updater_hom(),_1,_2),
					   boost::bind(multiplicative_phenotype_updater_het(),_1,_2),
					   1.);
  }
};

struct multiplicative_disease_effect_to_fitness
{
  typedef double result_type;
  template< typename iterator_type >
  inline double operator()(const iterator_type & g1, const iterator_type & g2,
			   const double & sd, const double & sd_s,gsl_rng * r) const
  {
    double genetic = multiplicative_phenotype()(g1,g2);
    double noise = (sd>0.)?gsl_ran_gaussian(r,sd):0.;
    //Subtract 1 so that phenotype has mean 1 and std_dev sd_s
    double fitness = exp( (-1. * pow(std::abs(genetic+noise)-1.,2.))/(2.*pow(sd_s,2)) );
    return ( fitness );
  }
};

#endif
