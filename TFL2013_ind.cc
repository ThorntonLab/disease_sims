/*
  Simulation from Thornton, Foran, and Long (2013) PLoS Genetics.

  Re-implemented using individual-based sampling routines.

  Rewritten to use the fwdpp library v >= 0.2.0
*/
#include <fwdpp/diploid.hh>
#include <utility>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

#include <boost/function.hpp>

#include <boost/unordered_set.hpp>
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/program_options.hpp>

#include <mutation_with_age.hpp>

#include <fcntl.h>

using namespace std;
using namespace boost::iostreams;
using namespace boost::program_options;
using namespace KTfwd;

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
    return site_dependent_fitness()(g1,g2,
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
    double noise = gsl_ran_gaussian(r,sd);
    //Subtract 1 so that phenotype has mean 1 and std_dev sd_s
    double fitness = exp( (-1. * pow(std::abs(genetic+noise)-1.,2.))/(2.*pow(sd_s,2)) );
    return ( fitness );
  }
};


//The mutation model
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

struct simparams
{
  unsigned N,N2,ngens_burnin,ngens_evolve,ngens_evolve_growth,replicate_no,seed;
  double mu_disease,mu_neutral,littler,s,sd,sd_s;
  bool dist_effects,multiplicative ;//= atoi(argv[argument++]);
  string indexfile, hapfile, phenofile, effectsfile ;
  simparams(void);
};

simparams::simparams(void) : N(20000),N2(20000),
			     ngens_burnin(0),
			     ngens_evolve(160000),
			     ngens_evolve_growth(0),
			     replicate_no(0),
			     seed(0),
			     mu_disease(0.000125),
			     mu_neutral(0.00125),
			     littler(0.00125),
			     s(0.1),
			     sd(0.075),
			     sd_s(1),
			     dist_effects(true),
			     multiplicative(false),
			     indexfile(string()),
			     hapfile(string()),
			     phenofile(string()),
			     effectsfile(string())
{
}

simparams parse_command_line(const int & argc,
			     char ** argv);



int main(int argc, char ** argv)
{
  simparams params = parse_command_line(argc,argv);

  //Determine growth rate under exponential growth model.
  double G = 0.;
  if( params.ngens_evolve_growth > 0 ) 
    { 
      G = exp( (log(double(params.N2)) - log(double(params.N)))/double(params.ngens_evolve_growth)); 
    }

  cerr << '#';
  for(int i=0;i<argc;++i)
    {
      cerr << argv[i] << ' ';
    }
  cerr << endl;

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,params.seed);

  lookup_table_type lookup;
  //the population begins with 1 gamete with no mutations
  glist gametes(1,gtype(2*params.N));
  mlist mutations;
  mvector fixations;      
  ftvector fixation_times;
  std::vector< std::pair< glist::iterator, glist::iterator> > diploids(params.N,
								       std::make_pair(gametes.begin(),
										      gametes.begin()));
  unsigned generation;
  unsigned ttl_gen = 0;
  double wbar=1;

  boost::function<double(void)> recmap = boost::bind(gsl_rng_uniform,r);

  for( generation = 0; generation < params.ngens_burnin; ++generation,++ttl_gen )
    {
      //Evolution w/no deleterious mutations and no selection.
      wbar = sample_diploid(r,
			    &gametes,
			    &diploids,
			    &mutations,
			    params.N,
			    params.mu_neutral,
			    boost::bind(mutation_model(),r,ttl_gen,params.s,0.,params.mu_neutral,&mutations,&lookup,params.dist_effects),
			    boost::bind(KTfwd::genetics101(),_1,_2,
					&gametes,
					params.littler,
					r,
					recmap),
			    boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
			    boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
			    boost::bind(KTfwd::no_selection(),_1,_2),
			    boost::bind(KTfwd::mutation_remover(),_1,0,2*params.N));
      KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,ttl_gen,2*params.N);
    }

  if (! params.multiplicative )
    {
      for( generation = 0; generation < params.ngens_evolve; ++generation,++ttl_gen )
	{
	  //Evolve under the disease model from TFL2013
	  wbar = sample_diploid(r,
				&gametes,
				&diploids,
				&mutations,
				params.N,
				params.mu_disease+params.mu_neutral,
				boost::bind(mutation_model(),r,ttl_gen,params.s,params.mu_disease,params.mu_neutral,&mutations,&lookup,params.dist_effects),
				boost::bind(KTfwd::genetics101(),_1,_2,
					    &gametes,
					    params.littler,
					    r,
					    recmap),
				boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
				boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
				boost::bind(disease_effect_to_fitness(),_1,_2,params.sd,params.sd_s,r),
				boost::bind(KTfwd::mutation_remover(),_1,0,2*params.N));
	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,ttl_gen,2*params.N);
	}
      //Exp. growth phase w/disease model
      for( generation = 0 ; generation < params.ngens_evolve_growth ; ++generation,++ttl_gen )
	{
	  unsigned N_next = round( params.N*pow(G,generation+1) );
	  wbar = sample_diploid(r,
				&gametes,
				&diploids,
				&mutations,
				params.N,
				N_next,
				params.mu_disease+params.mu_neutral,
				boost::bind(mutation_model(),r,ttl_gen,params.s,params.mu_disease,params.mu_neutral,&mutations,&lookup,params.dist_effects),
				boost::bind(KTfwd::genetics101(),_1,_2,
					    &gametes,
					    params.littler,
					    r,
					    recmap),
				boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
				boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
				boost::bind(disease_effect_to_fitness(),_1,_2,params.sd,params.sd_s,r),
				boost::bind(KTfwd::mutation_remover(),_1,0,2*params.N));
	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,ttl_gen,2*params.N);
	}
    }
  else //multiplicative phenotype model w/Gaussian stabilizing selection
    {
      for( generation = 0; generation < params.ngens_evolve; ++generation,++ttl_gen )
	{
	  wbar = sample_diploid(r,
				&gametes,
				&diploids,
				&mutations,
				params.N,
				params.mu_disease+params.mu_neutral,
				boost::bind(mutation_model(),r,ttl_gen,params.s,params.mu_disease,params.mu_neutral,&mutations,&lookup,params.dist_effects),
				boost::bind(KTfwd::genetics101(),_1,_2,
					    &gametes,
					    params.littler,
					    r,
					    recmap),
				boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
				boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
				boost::bind(multiplicative_disease_effect_to_fitness(),_1,_2,params.sd,params.sd_s,r),
				boost::bind(KTfwd::mutation_remover(),_1,0,2*params.N));
	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,ttl_gen,2*params.N);
	}
      for( generation = 0 ; generation < params.ngens_evolve_growth ; ++generation,++ttl_gen )
	{
	  unsigned N_next = round( params.N*pow(G,generation+1) );
	  wbar = sample_diploid(r,
				&gametes,
				&diploids,
				&mutations,
				params.N,
				N_next,
				params.mu_disease+params.mu_neutral,
				boost::bind(mutation_model(),r,ttl_gen,params.s,params.mu_disease,params.mu_neutral,&mutations,&lookup,params.dist_effects),
				boost::bind(KTfwd::genetics101(),_1,_2,
					    &gametes,
					    params.littler,
					    r,
					    recmap),
				boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
				boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
				boost::bind(multiplicative_disease_effect_to_fitness(),_1,_2,params.sd,params.sd_s,r),
				boost::bind(KTfwd::mutation_remover(),_1,0,2*params.N));
	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,ttl_gen,2*params.N);
	}
    }


  //Write out the population
  ostringstream popbuffer;
  write_binary_pop(&gametes,&mutations,&diploids,boost::bind(mwriter(),_1,_2),popbuffer);

  //Write out the phenotypes
  ostringstream phenobuffer;
  for( unsigned i = 0 ; i < diploids.size() ; ++i )
    {
      if (! params.multiplicative) //TFL disease model
	{
	  pair<double,double> pheno = disease_effect()(diploids[i].first,
						       diploids[i].second,
						       params.sd,
						       r);
	  double x = pheno.first;
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	  x = pheno.second;
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	}
      else //Multiplicative model
	{
	  double x = multiplicative_phenotype()(diploids[i].first,
						diploids[i].second);
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	  x = gsl_ran_gaussian(r,params.sd);
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	}
    }

  //Write out effects information for causative sites
  ostringstream effectstream;
  unsigned ncausative=0;
  for( typename mlist::const_iterator i = mutations.begin() ; i != mutations.end() ; ++i )
    {
      if ( ! i->neutral )
	{
	  ++ncausative;
	}
    }
  effectstream.write( reinterpret_cast<char *>(&ncausative),sizeof(unsigned) );
  for( typename mlist::const_iterator i = mutations.begin() ; i != mutations.end() ; ++i )
    {
      if(!i->neutral)
	{
	  if(! params.effectsfile.empty() )
	    {
	      //position of causative mutation, effect size, count, age
	      double pos = i->pos;
	      double s = i->s;
	      double count = double(i->n);
	      double age = double(ttl_gen - i->o + 1);
	      effectstream.write( reinterpret_cast<char *>(&pos), sizeof(double) );
	      effectstream.write( reinterpret_cast<char *>(&s), sizeof(double) );
	      effectstream.write( reinterpret_cast<char *>(&count), sizeof(double) );
	      effectstream.write( reinterpret_cast<char *>(&age), sizeof(double) );
	    }
	}
    }

  /*
    OK, now we lock the index file first.  Then, the rest of the output files.

    Write the relevant info to the index files, then the output files.

    Then, unlock output files first, releasing lock on index last.
  */
  //get file descriptors and grab locks in desired order
  struct flock index_flock, hapfile_flock, phenofile_flock, effects_flock;
  index_flock.l_type = F_WRLCK;/*Write lock*/
  index_flock.l_whence = SEEK_SET;
  index_flock.l_start = 0;
  index_flock.l_len = 0;/*Lock whole file*/

  //lock index file ASAP
  FILE * index_fh = fopen(params.indexfile.c_str(),"a");
  int index_fd = fileno(index_fh);
  if ( index_fd == -1 ) 
    { 
      std::cerr << "ERROR: could not open " << params.indexfile << '\n';
      exit(10);
    }
  if (fcntl(index_fd, F_SETLKW,&index_flock) == -1) 
    {
      std::cerr << "ERROR: could not obtain lock on " << params.indexfile << '\n';
      exit(10);
    }

  hapfile_flock.l_type = F_WRLCK;/*Write lock*/
  hapfile_flock.l_whence = SEEK_SET;
  hapfile_flock.l_start = 0;
  hapfile_flock.l_len = 0;/*Lock whole file*/

  phenofile_flock.l_type = F_WRLCK;/*Write lock*/
  phenofile_flock.l_whence = SEEK_SET;
  phenofile_flock.l_start = 0;
  phenofile_flock.l_len = 0;/*Lock whole file*/

  effects_flock.l_type = F_WRLCK;/*Write lock*/
  effects_flock.l_whence = SEEK_SET;
  effects_flock.l_start = 0;
  effects_flock.l_len = 0;/*Lock whole file*/

  FILE * haps_fh = fopen(params.hapfile.c_str(),"a");
  int hapfile_fd = fileno(haps_fh);

  if ( hapfile_fd == -1 ) 
    { 
      std::cerr << "ERROR: could not open " << params.hapfile << '\n';
      exit(10);
    }
  if (fcntl(index_fd, F_SETLKW,&index_flock) == -1) 
    {
      std::cerr << "ERROR: could not obtain lock on " << params.hapfile << '\n';
      exit(10);
    }

  FILE * pheno_fh = fopen(params.phenofile.c_str(),"a");
  int pheno_fd = fileno(pheno_fh);
  if ( pheno_fd == -1 ) 
    { 
      std::cerr << "ERROR: could not open " << params.phenofile << '\n';
      exit(10);
    }
  if (fcntl(pheno_fd, F_SETLKW,&index_flock) == -1) 
    {
      std::cerr << "ERROR: could not obtain lock on " << params.phenofile << '\n';
      exit(10);
    }

  FILE * effect_fh = NULL;
  int effect_fd = 0;
  effect_fh = fopen(params.effectsfile.c_str(),"a");
  effect_fd = fileno(effect_fh);
  if ( effect_fd == -1 ) 
    { 
      std::cerr << "ERROR: could not open " << params.effectsfile << '\n';
      exit(10);
    }
  if (fcntl(effect_fd, F_SETLKW,&index_flock) == -1) 
    {
      std::cerr << "ERROR: could not obtain lock on " << params.effectsfile << '\n';
      exit(10);
    }

  //Write the index data
  std::ostringstream indexstream;
  indexstream << params.replicate_no << ' ' << ftell(effect_fh) << ' '
	      << ftell(pheno_fh) << ' ' << ftell(haps_fh);
  fprintf(index_fh,"%s\n",indexstream.str().c_str());


  //write the haplotype data
  ::write(hapfile_fd,popbuffer.str().c_str(),popbuffer.str().size());
  //write the phenotype data
  ::write(pheno_fd, phenobuffer.str().c_str(), phenobuffer.str().size() );
  //write the effects
  ::write( effect_fd, effectstream.str().c_str(), effectstream.str().size() );
  //clear buffer

  //release the locks
  index_flock.l_type = F_UNLCK;
  hapfile_flock.l_type = F_UNLCK;
  phenofile_flock.l_type = F_UNLCK;
  effects_flock.l_type = F_UNLCK;

  if (fcntl(effect_fd, F_UNLCK,&effects_flock) == -1) 
    {
      std::cerr << "ERROR: could not release lock on " << params.effectsfile << '\n';
      exit(10);
    }
  fflush(effect_fh);
  fclose(effect_fh);

  if (fcntl(pheno_fd, F_UNLCK,&phenofile_flock) == -1) 
    {
      std::cerr << "ERROR: could not release lock on " << params.phenofile << '\n';
      exit(10);
    }
  fflush( pheno_fh );
  fclose( pheno_fh);

  if (fcntl(hapfile_fd, F_UNLCK,&hapfile_flock) == -1) 
    {
      std::cerr << "ERROR: could not releaselock on " <<  params.hapfile << '\n';
      exit(10);
    }
  fflush( haps_fh );
  fclose(haps_fh);

  if (fcntl(index_fd, F_UNLCK,&index_flock) == -1) 
    {
      std::cerr << "ERROR: could not release lock on " << params.indexfile << '\n';
      exit(10);
    }
  fflush( index_fh );
  fclose(index_fh);
}

simparams parse_command_line(const int & argc,
			     char ** argv)
{
  simparams rv;
  options_description desc("Simulate the genetic model from Thornton, Foran, and Long (2013) PLoS Genetics 9(2): e1003258.\nUsage: TFL2013 -h to see help");
  desc.add_options()
    ("help,h", "Produce help message")
    ("popsize1,1", value<unsigned>(&rv.N)->default_value(20000), "Diploid population number")
    ("popsize2,2", value<unsigned>(&rv.N2)->default_value(20000), "Population size to grow towards")
    ("burnin",value<unsigned>(&rv.ngens_burnin)->default_value(0), "No. generations to evolve with neutral mutations and no selection")
    ("generations,g",value<unsigned>(&rv.ngens_evolve)->default_value(160000), "No. generations to evolve with neutral mutations, deleterious mutations, and selection")
    ("generation-growth,G",value<unsigned>(&rv.ngens_evolve_growth)->default_value(0), "No. generation to grow exponentially from population size popsize1 to popsize2. Default is 0, meaning growth does not happen")
    ("replicate,R",value<unsigned>(&rv.replicate_no)->default_value(0),"Label of this replicate.  Must be >= 0")
    ("neutral,n",value<double>(&rv.mu_neutral)->default_value(0.00125),"Neutral mutation rate (per region per generation")
    ("causative,c",value<double>(&rv.mu_neutral)->default_value(0.000125),"Mutation rate to causative mutations(per region per generation")
    ("recrate,r",value<double>(&rv.littler)->default_value(0.00125),"Recombination rate (per diploid per region per generation")
    ("esize,e",value<double>(&rv.s)->default_value(0.1),"Effect size of causative mutation.  Mean of exponential dist by default.  Constant effect size if dist 0 or -d 0 is used")
    ("noise",value<double>(&rv.sd)->default_value(0.075),"Std. deviation in Gaussian noise to be added to phenotype")
    ("sigma",value<double>(&rv.sd_s)->default_value(1.0),"Std. deviation in Gaussian fitness function")
    ("dist,d",value<bool>(&rv.dist_effects)->default_value(true),"If true, model distribution of effect sizes.  Otherwise, constant effect size")
    ("multiplicative,m",value<bool>(&rv.multiplicative)->default_value(false),"If true, use multiplicative phenotype model.  Default is Thornton, Foran & Long (2013) model")
    ("indexfile,i",value<string>(&rv.indexfile)->default_value(string()),"Name of index file")
    ("popfile,p",value<string>(&rv.hapfile)->default_value(string()),"Name of output file for population")
    ("phenotypes,P",value<string>(&rv.phenofile)->default_value(string()),"Name of output file for phenotypes")
    ("effectsfile,E",value<string>(&rv.effectsfile)->default_value(string()),"Name of output file for effect sizes of causative mutations")
    ("seed,S",value<unsigned>(&rv.seed)->default_value(0),"Random number seed (unsigned integer)")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if(argc == 1 || vm.count("help"))
    {
      cerr << desc << '\n';
      exit(0);
    }

  if( rv.indexfile.empty() || rv.hapfile.empty() || rv.phenofile.empty() )
    {
      if (rv.indexfile.empty())
	{
	  cerr << "Error: indef file name required.  Use -h to see options\n";
	}

      if( rv.hapfile.empty() )
	{
	  cerr << "Error: population output file name required.  Use -h to see options\n";
	}
	
      if( rv.phenofile.empty() )
	{
	  cerr << "Error: phenotypes output file name required.  Use -h to see options\n";
	}
      exit(10);
    }

  return rv;
}
