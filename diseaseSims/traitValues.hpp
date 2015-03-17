#ifndef __TRAIT_VALUES_HPP__
#define __TRAIT_VALUES_HPP__

#include <fwdpp/fitness_models.hpp>

/*
  Functions for TFL2013 and additive trait calculators.
  The simulation uses calculators that return pairs of doubles.
  The pair is the G,E components of trait value.
  Here, we just need the G, so we do this as a quick fix:
*/
struct TFL2013g
{
  using return_type = double;
  inline double operator()( const glist::const_iterator & g1,
			    const glist::const_iterator & g2) const
  {
    double e1 = std::accumulate( g1->smutations.begin(),
				 g1->smutations.end(),
				 0.,
				 [](const double & a,
				    const gtype::mutation_list_type_iterator & b) { return a + b->s; } );
    double e2 = std::accumulate( g2->smutations.begin(),
				 g2->smutations.end(),
				 0.,
				 [](const double & a,
				    const gtype::mutation_list_type_iterator & b) { return a + b->s; } );
    return (sqrt(e1*e2));
  }
};

struct additiveg
{
  using return_type = double;
  inline double operator()( const glist::const_iterator & g1,
			    const glist::const_iterator & g2) const
  {
    return std::accumulate( g1->smutations.begin(),
			    g1->smutations.end(),
			    0.,
			    [](const double & a,
			       const gtype::mutation_list_type_iterator & b) { return a + b->s; } )
      + std::accumulate( g2->smutations.begin(),
			 g2->smutations.end(),
			 0.,
			 [](const double & a,
			    const gtype::mutation_list_type_iterator & b) { return a + b->s; } );
  }
};

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

#endif
