/*
  Simulation from Thornton, Foran, and Long (2013) PLoS Genetics.

  Rewritten to use the fwdpp library v >= 0.2.0
*/
#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <utility>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

#ifndef USE_STANDARD_CONTAINERS
#include <boost/unordered_set.hpp>
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#else
#include <vector>
#include <list>
#include <unordered_set>
#include <functional>
#endif
#include <iostream>
#include <sstream>
using namespace std;
using namespace std::placeholders;
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
typedef boost::container::vector<gtype,gam_allocator > gvector;
typedef boost::container::vector<mtype> mvector;
typedef boost::container::vector<unsigned> ftvector;
typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;
#else
typedef std::list<mtype > mlist;
typedef gamete_base<mtype, mlist> gtype;
typedef vector<gtype> gvector;
typedef vector<mtype> mvector;
typedef vector<unsigned> ftvector;
typedef std::unordered_set<double,std::hash<double>,KTfwd::equal_eps > lookup_table_type;
#endif



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
  const double littler = atof(argv[argument++])/2.; //convert per-diploid into per-gamete
  const double s = atof(argv[argument++]);
  const bool dist_effects = atoi(argv[argument++]);
  const double sd = atof(argv[argument++]);//sigma for "E"
  const double sd_s = atof(argv[argument++]);//sigma for selection
  const unsigned ngens_burnin = atoi(argv[argument++]);
  const unsigned ngens_evolve = atoi(argv[argument++]);
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

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);

  while(nreps--)
    {
      lookup_table_type lookup;
      //the population begins with 1 gamete with no mutations
      gvector gametes(1,gtype(2*N));
      mlist mutations;
      mvector fixations;      
      ftvector fixation_times;

      unsigned generation;
      unsigned ttl_gen = 0;
      double wbar=1;
      for( generation = 0; generation < ngens_burnin; ++generation,++ttl_gen )
	{
	  //a generation is drift, then mutation, recombination
	  wbar = sample_diploid(r,&gametes,2*N,std::bind(disease_effect_to_fitness(),std::placeholders::_1,std::placeholders::_2,sd,sd_s,r),
				std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*N));       
	  remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,2*N);
	  assert(check_sum(gametes,2*N));
	  unsigned nmuts= mutate(r,&gametes,&mutations,mu_disease+mu_neutral,std::bind(mutation_model(),r,ttl_gen,s,mu_disease,mu_neutral,&mutations,&lookup,dist_effects),
				 push_back_gamete(),std::bind(insert_mutation_at_end<mtype,mlist >,std::placeholders::_1,std::placeholders::_2));
	  assert(check_sum(gametes,2*N));
	  unsigned nrec = recombine(r, &gametes, 2*N, littler, std::bind(gsl_rng_uniform,r));
	  assert(check_sum(gametes,2*N));
	}
      unsigned nfix = fixations.size();
      for( generation = 0; generation < ngens_evolve; ++generation,++ttl_gen )
	{
	  //a generation is drift, then mutation, recombination
	  double wbar = sample_diploid(r,&gametes,2*N,std::bind(disease_effect_to_fitness(),std::placeholders::_1,std::placeholders::_2,sd,sd_s,r),
				       std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*N));       
	  remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,2*N);
	  assert(check_sum(gametes,2*N));
	  unsigned nmuts= mutate(r,&gametes,&mutations,mu_disease+mu_neutral,std::bind(mutation_model(),r,ttl_gen,s,mu_disease,mu_neutral,&mutations,&lookup,dist_effects),
				 push_back_gamete(),std::bind(insert_mutation_at_end<mtype,mlist >,std::placeholders::_1,std::placeholders::_2));
	  assert(check_sum(gametes,2*N));
	  unsigned nrec = recombine(r, &gametes, 2*N, littler, std::bind(gsl_rng_uniform,r));
	  assert(check_sum(gametes,2*N));
	}
      remove_fixed_lost(&mutations,&fixations,&fixation_times,generation,2*N);

      //so now, make 2N diploids from the population as it stands, and assign a phenotype value to each individual
      vector< pair< gvector::iterator,gvector::iterator > > diploids;
      vector< unsigned > gamcounts;
      vector< double > phenotypes;
      unsigned gamsum=0;
      for(unsigned i=0;i<gametes.size();++i)
	{
	  gamcounts.push_back(gametes[i].n);
	  gamsum += gametes[i].n;
	}

      ostringstream buffer;
      //generate the 2N diploids
      unsigned gam1,gam2,j,k,l;
      for(unsigned i=0;i<N;)
	{
	  gam1 = unsigned(gsl_ran_flat(r,0,gamsum))+1;
	  gam2 = unsigned(gsl_ran_flat(r,0,gamsum))+1;
	  j=0;
	  l=0;
	  for( ;l < gamcounts.size() ; ++l)
	    {
	      j += gamcounts[l];
	      if( j >= gam1 ) { j=l;break; }
	    }
	  k=0;
	  l=0;
	  for( ;l < gamcounts.size() ; ++l)
	    {
	      k += gamcounts[l];
	      if( k >= gam2 ) { k=l;break; }
	    }
 	  if( ( k==j && gamcounts[k]>1 ) ||
 	      ( k!=j && gamcounts[k]>0 && gamcounts[j]>0 ) )
 	    {
	      diploids.push_back( std::make_pair( gametes.begin()+j, gametes.begin()+k ) );
	      pair<double,double> effect = disease_effect()(gametes.begin()+j,gametes.begin()+k,sd,r);
	      buffer << effect.first << '\t' << effect.second << '\n';
	      gamcounts[j]--;
	      gamcounts[k]--;
	      gamsum -= 2;
	      ++i;
	    }
	}
      assert( accumulate(gamcounts.begin(),gamcounts.end(),0u) == 0 );
      ofstream effectfile(ofn2);
      effectfile << buffer.str();
      effectfile.close();
      buffer.str(string());
      assert( gamsum == 0 );

      //print out num deleterious mutations vs effect
      vector< pair<double,string> > neutral,selected;
      SimData msneut,mssel;
      std::function<bool(const std::pair<double,std::string> &, const double &)> sitefinder = [](const std::pair<double,std::string> & site,
												 const double & d ) 
	{
	  return std::fabs(site.first-d) <= std::numeric_limits<double>::epsilon();
	};
      for(unsigned i=0;i<diploids.size();++i)
	{
	  //neutral mutations
	  for( unsigned mut = 0 ; mut < diploids[i].first->mutations.size() ; ++mut )
	    {
	      vector< pair<double,string> >::iterator itr = find_if(neutral.begin(),neutral.end(),std::bind(sitefinder, std::placeholders::_1,diploids[i].first->mutations[mut]->pos));
	      if( itr == neutral.end() )
		{
		  neutral.push_back( std::make_pair(diploids[i].first->mutations[mut]->pos,std::string(2*N,'0')) );
		  neutral[neutral.size()-1].second[2*i] = '1';
		}
	      else
		{
		  itr->second[2*i] = '1';
		}
	    }
	  for( unsigned mut = 0 ; mut < diploids[i].second->mutations.size() ; ++mut )
	    {
	      vector< pair<double,string> >::iterator itr = find_if(neutral.begin(),neutral.end(),std::bind(sitefinder, std::placeholders::_1, diploids[i].second->mutations[mut]->pos));
	      if( itr == neutral.end() )
		{
		  neutral.push_back( std::make_pair(diploids[i].second->mutations[mut]->pos,std::string(2*N,'0')) );
		  neutral[neutral.size()-1].second[2*i+1] = '1';
		}
	      else
		{
		  itr->second[2*i+1] = '1';
		}
	    }

	  //selected
	  for( unsigned mut = 0 ; mut < diploids[i].first->smutations.size() ; ++mut )
	    {
	      vector< pair<double,string> >::iterator itr = find_if(selected.begin(),selected.end(),std::bind(sitefinder,std::placeholders::_1, diploids[i].first->smutations[mut]->pos));
	      if( itr == selected.end() )
		{
		  selected.push_back( std::make_pair(diploids[i].first->smutations[mut]->pos,std::string(2*N,'0')) );
		  selected[selected.size()-1].second[2*i] = '1';
		}
	      else
		{
		  itr->second[2*i] = '1';
		}
	    }
	  for( unsigned mut = 0 ; mut < diploids[i].second->smutations.size() ; ++mut )
	    {
	      vector< pair<double,string> >::iterator itr = find_if(selected.begin(),selected.end(),std::bind(sitefinder,std::placeholders::_1, diploids[i].second->smutations[mut]->pos));
	      if( itr == selected.end() )
		{
		  selected.push_back( make_pair(diploids[i].second->smutations[mut]->pos,std::string(2*N,'0')) );
		  selected[selected.size()-1].second[2*i+1] = '1';
		}
	      else
		{
		  itr->second[2*i+1] = '1';
		}
	    }
	}
      std::sort(neutral.begin(),neutral.end(),
		[](std::pair<double,std::string> lhs,
		   std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
      std::sort(selected.begin(),selected.end(),
		[](std::pair<double,std::string> lhs,
		   std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
      msneut.assign(neutral.begin(),neutral.end());
      mssel.assign(selected.begin(),selected.end());

      filtering_ostream msfile;
      msfile.push(gzip_compressor());
      msfile.push(file_sink(ofn,ios_base::out|ios_base::binary));
      msfile << msneut << '\n'<<mssel<<'\n';

      buffer << "//\n";
      for(mlist::iterator i = mutations.begin() ; i != mutations.end() ; ++i )
	{
	  if (!i->neutral)
	    {
	      buffer << i->pos << '\t' << i->s << '\n';
	    }
	}
      filtering_ostream mutationfile;
      mutationfile.push(gzip_compressor());
      mutationfile.push(file_sink(ofn3,ios_base::out|ios_base::binary));
      mutationfile << buffer.str();
    }
}
