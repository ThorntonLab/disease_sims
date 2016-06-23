#ifndef __TRAIT_VALUES_HPP__
#define __TRAIT_VALUES_HPP__

#include <numeric>
#include <fwdpp/fitness_models.hpp>
#include <diseaseSims/mutation_with_age.hpp>
/*
  Geometric mean of maternal/paternal haplotypes.
*/
struct TFL2013g
{
  using return_type = double;
  inline double operator()( const poptype::gamete_t & g1,
			    const poptype::gamete_t & g2,
			    const poptype::mcont_t & mutations) const
  {
    double e1 = std::accumulate( g1.smutations.begin(),
				 g1.smutations.end(),
				 0.,
				 [&mutations](const double & a,
					      const std::size_t b) { return a + mutations[b].s; } );
    double e2 = std::accumulate( g2.smutations.begin(),
				 g2.smutations.end(),
				 0.,
				 [&mutations](const double & a,
					      const std::size_t b) { return a + mutations[b].s; } );
    return (sqrt(e1*e2));
  }
};

//Additive
struct additiveg
{
  using return_type = double;
  inline double operator()( const poptype::gamete_t & g1,
			    const poptype::gamete_t & g2,
			    const poptype::mcont_t & mutations) const
  {
    return std::accumulate( g1.smutations.begin(),
			    g1.smutations.end(),
			    0.,
			    [&mutations](const double & a,
					 const std::size_t b) { return a + mutations[b].s; } )
      + std::accumulate( g2.smutations.begin(),
			 g2.smutations.end(),
			 0.,
			 [&mutations](const double & a,
				      const std::size_t b) { return a + mutations[b].s; } );
  }
};

//multiplicative
struct multiplicative_phenotype
{
  typedef double result_type;
  template< typename gamete_type,
	    typename mcont_t>
  inline double operator()(const gamete_type & g1, const gamete_type & g2,const mcont_t & mutations) const
  {
    using __mtype =  typename mcont_t::value_type;
    return KTfwd::site_dependent_fitness()(g1,g2,mutations,
					   [](double & fitness,const __mtype & mut)
					   {
					     fitness *= ( std::pow(1. + mut.s,2.) ); 
					   },
					   [](double & fitness,const __mtype & mut)
					   {
					     fitness *= ( 1. + mut.s ); 
					   },
					   1.);
  }
};

//Popgen-like model for trait values
struct popgen_phenotype_updater_hom
{
  typedef void result_type;
  template<typename mutation_type>
  inline void operator()(double & fitness, const mutation_type & m1) const
  {
    /*
      Bug fix March 4, 2016
      Old code:     fitness *= (1. + m1->s);
      It was a bug b/c the trait values
      didn't scale as they would for homozygote
      in Risch model.
    */
    fitness *= (1. + 2.*m1.s);
  }
};

struct popgen_phenotype_updater_het
{
  typedef void result_type;
  template<typename mutation_type>
  inline void operator()(double & fitness, const mutation_type & m1,const double & dominance) const
  {
    fitness *= ( 1. + dominance*m1.s );
  }
};

struct popgen_phenotype
{
  typedef double result_type;
  template< typename gamete_type,typename mcont_t>
  inline double operator()(const gamete_type & g1, const gamete_type & g2,const mcont_t & mutations,
			   const double & dominance) const
  {
    return KTfwd::site_dependent_fitness()(g1,g2,mutations,
					   std::bind(popgen_phenotype_updater_hom(),std::placeholders::_1,std::placeholders::_2),
					   std::bind(popgen_phenotype_updater_het(),std::placeholders::_1,std::placeholders::_2,dominance),
					   1.);
  }
};

#endif
