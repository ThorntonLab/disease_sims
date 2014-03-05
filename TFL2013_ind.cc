/*
  Simulation from Thornton, Foran, and Long (2013) PLoS Genetics.

  Re-implemented using individual-based sampling routines.

  Rewritten to use the fwdpp library v >= 0.2.0
*/
#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <utility>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

#include <boost/function.hpp>

#include <boost/unordered_set.hpp>
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>

using namespace std;
using namespace boost::iostreams;
using Sequence::SimData;
using namespace KTfwd;

struct mutation_with_age : public mutation_base
{
  mutable unsigned o;
  double s;
  char label; //'A' = "amino acid", 'S' = "synonymous"  Only used here in that A = causative, S = neutral.
  mutation_with_age( const double & position, const double & sel_coeff,
		     const unsigned & count,const unsigned & origin, const char & ch,
		     const bool & n=true) 
    : mutation_base(position,count,n),o(origin),s(sel_coeff),label(ch)
  {
  }
  bool operator==(const mutation_with_age & rhs) const
  {
    return( fabs(this->pos-rhs.pos) <= std::numeric_limits<double>::epsilon() &&
	    this->s == rhs.s );
  }	
};

typedef mutation_with_age mtype;
//boost containers
#ifndef USE_STANDARD_CONTAINERS
typedef boost::pool_allocator<mtype> mut_allocator;
typedef boost::container::list<mtype,mut_allocator > mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef boost::pool_allocator<gtype> gam_allocator;
typedef boost::container::list<gtype,gam_allocator > glist;
typedef boost::container::vector<mtype> mvector;
typedef boost::container::vector<unsigned> ftvector;
#else
typedef std::list<mtype > mlist;
typedef gamete_base<mtype, mlist> gtype;
typedef std::list<gtype> glist;
typedef vector<mtype> mvector;
typedef vector<unsigned> ftvector;
#endif

typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;

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
    return make_pair(effect, gsl_ran_gaussian(r,sd));
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
    pair<double,double> effect = disease_effect()(g1,g2,sd,r);
    double fitness = exp( (-1. * pow(effect.first+effect.second,2.))/(2.*pow(sd_s,2)) );
    return ( fitness );
  }
};

struct mutation_model
{
  typedef mtype result_type;
  inline result_type operator()( gsl_rng * r, const unsigned int & ttl_generations,
				 const double & s, const double & ud, const double & un, const  mlist * mutations,
				 lookup_table_type * lookup,
				 const bool dist_effects = false) const
  {
    double pos = gsl_rng_uniform(r);
    while(lookup->find(pos) != lookup->end())
      {
	pos = gsl_rng_uniform(r);
      }
    lookup->insert(pos);
    if( gsl_rng_uniform(r) <= ud/(ud+un) )
      {
	if( ! dist_effects )
	  {
	    return mtype(pos,s,1,ttl_generations,'A',false);
	  }
	else
	  {
	    return mtype(pos,gsl_ran_exponential(r,s),1,ttl_generations,'A',false);
	  }
      }
    return mtype(pos,0.,1,ttl_generations,'S',true);
  }
};

int main(int argc, char ** argv)
{
  int argument=1;
  const unsigned N = atoi(argv[argument++]);
  const double mu_disease = atof(argv[argument++]);
  const double mu_neutral = atof(argv[argument++]);
  const double littler = atof(argv[argument++]);
  const double s = atof(argv[argument++]);
  const bool dist_effects = atoi(argv[argument++]);
  const double sd = atof(argv[argument++]);//sigma for "E"
  const double sd_s = atof(argv[argument++]);//sigma for selection
  const unsigned ngens_burnin = atoi(argv[argument++]);
  const unsigned ngens_evolve = atoi(argv[argument++]);
  const unsigned N2 = atoi(argv[argument++]);
  const unsigned ngens_evolve_growth = atoi(argv[argument++]);
  const unsigned samplesize1 = atoi(argv[argument++]);
  const char * ofn = argv[argument++];
  const char * ofn2 = argv[argument++];
  const char * ofn3 = NULL;
  if( dist_effects )
    {
      ofn3 = argv[argument++];
    }
  int nreps=atoi(argv[argument++]);
  const unsigned seed = atoi(argv[argument++]);

  //Determine growth rate under exponential growth model.
  double G = 0.;
  if( ngens_evolve_growth > 0 ) 
    { 
      G = exp( (log(double(N2)) - log(double(N)))/double(ngens_evolve_growth)); 
    }

  cerr << '#';
  for(unsigned i=0;i<argc;++i)
    {
      cerr << argv[i] << ' ';
    }
  cerr << endl;

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);

  while(nreps--)
    {
      lookup_table_type lookup;
      //the population begins with 1 gamete with no mutations
      glist gametes(1,gtype(2*N));
      mlist mutations;
      mvector fixations;      
      ftvector fixation_times;
      std::vector< std::pair< glist::iterator, glist::iterator> > diploids(N,
									   std::make_pair(gametes.begin(),
											  gametes.begin()));
      unsigned generation;
      unsigned ttl_gen = 0;
      double wbar=1;

      boost::function<double(void)> recmap = boost::bind(gsl_rng_uniform,r);

      for( generation = 0; generation < ngens_burnin; ++generation,++ttl_gen )
	{
	  //Evolution w/no deleterious mutations and no selection.
	  wbar = sample_diploid(r,
				&gametes,
				&diploids,
				&mutations,
				N,
				mu_neutral,
				boost::bind(mutation_model(),r,ttl_gen,s,0.,mu_neutral,&mutations,&lookup,dist_effects),
				boost::bind(KTfwd::genetics101(),_1,_2,
					    &gametes,
					    littler,
					    r,
					    recmap),
				boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
				boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
				boost::bind(KTfwd::no_selection(),_1,_2),
				boost::bind(KTfwd::mutation_remover(),_1,0,2*N));
      	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,ttl_gen,2*N);
	}

      for( generation = 0; generation < ngens_evolve; ++generation,++ttl_gen )
	{
	  //Evolve under the disease model
	  wbar = sample_diploid(r,
				&gametes,
				&diploids,
				&mutations,
				N,
				mu_disease+mu_neutral,
				boost::bind(mutation_model(),r,ttl_gen,s,mu_disease,mu_neutral,&mutations,&lookup,dist_effects),
				boost::bind(KTfwd::genetics101(),_1,_2,
					    &gametes,
					    littler,
					    r,
					    recmap),
				boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
				boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
				boost::bind(disease_effect_to_fitness(),_1,_2,sd,sd_s,r),
				boost::bind(KTfwd::mutation_remover(),_1,0,2*N));
      	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,ttl_gen,2*N);
	}
      //Exp. growth phase w/disease model
      for( generation = 0 ; generation < ngens_evolve_growth ; ++generation,++ttl_gen )
	{
	  unsigned N_next = round( N*pow(G,generation+1) );
	  wbar = sample_diploid(r,
				&gametes,
				&diploids,
				&mutations,
				N,
				N_next,
				mu_disease+mu_neutral,
				boost::bind(mutation_model(),r,ttl_gen,s,mu_disease,mu_neutral,&mutations,&lookup,dist_effects),
				boost::bind(KTfwd::genetics101(),_1,_2,
					    &gametes,
					    littler,
					    r,
					    recmap),
				boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
				boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
				boost::bind(disease_effect_to_fitness(),_1,_2,sd,sd_s,r),
				boost::bind(KTfwd::mutation_remover(),_1,0,2*N));
      	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,ttl_gen,2*N);
	}
    }
}
