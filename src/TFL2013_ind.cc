/*
  Simulation from Thornton, Foran, and Long (2013) PLoS Genetics.

  Re-implemented using individual-based sampling routines.

  Rewritten to use the fwdpp library v >= 0.2.4
*/
#include <fwdpp/diploid.hh>
#include <utility>
#include <iostream>

#include <boost/program_options.hpp>

#include <fcntl.h>

#include <zlib.h>

#include <mutation_with_age.hpp>
#include <TFL_fitness_models.hpp>

using namespace std;
using namespace boost::program_options;
using namespace KTfwd;

//The mutation model
struct mutation_model
{
  typedef TFLmtype result_type;
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
	    return TFLmtype(pos,s,1,ttl_generations,'A',false);
	  }
	else
	  {
	    return TFLmtype(pos,gsl_ran_exponential(r,s),1,ttl_generations,'A',false);
	  }
      }
    return TFLmtype(pos,0.,1,ttl_generations,'S',true);
  }
};

//fitness models
enum MODEL { GENE_RECESSIVE, GENE_ADDITIVE, MULTIPLICATIVE, POPGEN };

struct simparams
{
  unsigned N,N2,ngens_burnin,ngens_evolve,ngens_evolve_growth,replicate_no,seed;
  double mu_disease,mu_neutral,littler,s,sd,sd_s,optimum,dominance;
  bool dist_effects,gzoutput;
  MODEL model;
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
			     optimum(0.),
			     dominance(0.),
			     dist_effects(true),
			     gzoutput(false),
			     model( GENE_RECESSIVE ),
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

  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r);

  for( generation = 0; generation < params.ngens_burnin; ++generation,++ttl_gen )
    {
      //Evolution w/no deleterious mutations and no selection.
      wbar = sample_diploid(r,
			    &gametes,
			    &diploids,
			    &mutations,
			    params.N,
			    params.mu_neutral,
			    std::bind(mutation_model(),r,ttl_gen,params.s,0.,params.mu_neutral,&mutations,&lookup,params.dist_effects),
			    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
					&gametes,
					params.littler,
					r,
					recmap),
			    std::bind(KTfwd::insert_at_end<TFLmtype,mlist>,std::placeholders::_1,std::placeholders::_2),
			    std::bind(KTfwd::insert_at_end<gtype,glist>,std::placeholders::_1,std::placeholders::_2),
			    std::bind(KTfwd::no_selection(),std::placeholders::_1,std::placeholders::_2),
			    std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*params.N));
      KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,ttl_gen,2*params.N);
    }
  unsigned N_current = params.N;

  //Fitness model for phase w/selection.  The default is the recessive model of TFL 2013
  std::function<double(const glist::const_iterator &,
		       const glist::const_iterator &)> dipfit = std::bind(TFL2013_recessive(),std::placeholders::_1,std::placeholders::_2,params.sd,params.sd_s,0.,r);

  if( params.model == GENE_ADDITIVE )
    {
      dipfit = std::bind(TFL2013_additive(),std::placeholders::_1,std::placeholders::_2,params.sd,params.sd_s,0.,r);
    }
  else if( params.model == MULTIPLICATIVE )
    {
      dipfit = std::bind(multiplicative_disease_effect_to_fitness(),std::placeholders::_1,std::placeholders::_2,
			 params.sd,params.sd_s,params.optimum,r);
    }
  else if ( params.model == POPGEN )
    {
      dipfit = std::bind(popgen_disease_effect_to_fitness(),std::placeholders::_1,std::placeholders::_2,params.dominance,
			   params.sd,params.sd_s,params.optimum,r);
    }

  for( generation = 0; generation < params.ngens_evolve; ++generation,++ttl_gen )
    {
      //Evolve under the disease model from TFL2013
      wbar = sample_diploid(r,
			    &gametes,
			    &diploids,
			    &mutations,
			    params.N,
			    params.mu_disease+params.mu_neutral,
			    std::bind(mutation_model(),r,ttl_gen,params.s,params.mu_disease,params.mu_neutral,&mutations,&lookup,params.dist_effects),
			    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
					&gametes,
					params.littler,
					r,
					recmap),
			    std::bind(KTfwd::insert_at_end<TFLmtype,mlist>,std::placeholders::_1,std::placeholders::_2),
			    std::bind(KTfwd::insert_at_end<gtype,glist>,std::placeholders::_1,std::placeholders::_2),
			    dipfit,
			    std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*params.N));
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
			    N_current,
			    N_next,
			    params.mu_disease+params.mu_neutral,
			    std::bind(mutation_model(),r,ttl_gen,params.s,params.mu_disease,params.mu_neutral,&mutations,&lookup,params.dist_effects),
			    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
					&gametes,
					params.littler,
					r,
					recmap),
			    std::bind(KTfwd::insert_at_end<TFLmtype,mlist>,std::placeholders::_1,std::placeholders::_2),
			    std::bind(KTfwd::insert_at_end<gtype,glist>,std::placeholders::_1,std::placeholders::_2),
			    dipfit,
			    std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*N_next));
      KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,ttl_gen,2*N_next);
      //update N
      N_current = N_next;
    }

  //Write out the population
  ostringstream popbuffer;
  write_binary_pop(&gametes,&mutations,&diploids,std::bind(mwriter(),std::placeholders::_1,std::placeholders::_2),popbuffer);

  //Write out the phenotypes
  ostringstream phenobuffer;
  unsigned nphenos = diploids.size();
  phenobuffer.write( reinterpret_cast<char *>(&nphenos), sizeof(unsigned) );
  for( unsigned i = 0 ; i < diploids.size() ; ++i )
    {
      if (params.model == GENE_RECESSIVE) 
	{
	  pair<double,double> pheno = TFL2013_recessive_disease_effect()(diploids[i].first,
									 diploids[i].second,
									 params.sd,
									 r);
	  double x = pheno.first;
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	  x = pheno.second;
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	}
      if (params.model == GENE_ADDITIVE) 
	{
	  pair<double,double> pheno = TFL2013_additive_disease_effect()(diploids[i].first,
									diploids[i].second,
									params.sd,
									r);
	  double x = pheno.first;
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	  x = pheno.second;
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	}
      else if (params.model == MULTIPLICATIVE) 
	{
	  double x = multiplicative_phenotype()(diploids[i].first,
						diploids[i].second);
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	  x = (params.sd > 0.) ? gsl_ran_gaussian(r,params.sd) : 0.;
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	}
      else if (params.model == POPGEN)
	{
	  double x = popgen_phenotype()(diploids[i].first,
					diploids[i].second,params.dominance);
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	  x = (params.sd > 0.) ? gsl_ran_gaussian(r,params.sd) : 0.;
	  phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	}
    }

  //Write out effects information for causative sites
  ostringstream effectstream;
  unsigned nmuts = mutations.size();
  /*
    unsigned ncausative=0;
    for( typename mlist::const_iterator i = mutations.begin() ; i != mutations.end() ; ++i )
    {
    if ( ! i->neutral )
    {
    ++ncausative;
    }
    }
  */
  //effectstream.write( reinterpret_cast<char *>(&ncausative),sizeof(unsigned) );
  effectstream.write( reinterpret_cast<char *>(&nmuts),sizeof(unsigned) );
  for( typename mlist::const_iterator i = mutations.begin() ; i != mutations.end() ; ++i )
    {
      //if(!i->neutral)
      //{
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
      //}
    }
  /*
    OK, now we lock the index file first.  Then, the rest of the output files.

    Write the relevant info to the index files, then the output files.

    Then, unlock output files first, releasing lock on index last.
  */
  //get file descriptors and grab locks in desired order
  struct flock index_flock;//, hapfile_flock, phenofile_flock, effects_flock;
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

  if( params.gzoutput ) 
    {
      gzFile gzout = gzopen(params.hapfile.c_str(),"a");
      int hapswritten = gzwrite(gzout,popbuffer.str().c_str(),popbuffer.str().size());
      if ( ! hapswritten && !popbuffer.str().empty() )
	{
	  cerr << "Error writing population to " 
	       << params.hapfile << '\n';
	  exit(10);
	}
      //gzflush(gzout,Z_FINISH);
      gzclose(gzout);

      gzout = gzopen(params.phenofile.c_str(),"a");
      int phenowritten = gzwrite(gzout,phenobuffer.str().c_str(),phenobuffer.str().size());
      if ( ! phenowritten && !phenobuffer.str().empty() )
	{
	  cerr << "Error writing population to " 
	       << params.phenofile << '\n';
	  exit(10);
	}
      //gzflush(gzout,Z_FINISH);
      gzclose(gzout);

      gzout = gzopen(params.effectsfile.c_str(),"a");
      int effectwritten = gzwrite(gzout,effectstream.str().c_str(),effectstream.str().size());
      if ( ! effectwritten && !effectstream.str().empty() )
	{
	  cerr << "Error writing population to " 
	       << params.effectsfile << '\n';
	  exit(10);
	}
      //gzflush(gzout,Z_FINISH);
      gzclose(gzout);

      //Now we can write to the index file
      std::ostringstream indexstream;
      indexstream << params.replicate_no << ' ' << effectwritten << ' '
		  << phenowritten << ' ' << hapswritten;
      fprintf(index_fh,"%s\n",indexstream.str().c_str());
    }
  else
    {
      FILE * haps_fh = fopen(params.hapfile.c_str(),"a");
      int hapfile_fd = fileno(haps_fh);
      if ( hapfile_fd == -1 ) 
	{ 
	  std::cerr << "ERROR: could not open " << params.hapfile << '\n';
	  exit(10);
	}
      FILE * pheno_fh = fopen(params.phenofile.c_str(),"a");
      int pheno_fd = fileno(pheno_fh);
      if ( pheno_fd == -1 ) 
	{ 
	  std::cerr << "ERROR: could not open " << params.phenofile << '\n';
	  exit(10);
	}
      FILE * effect_fh = fopen(params.effectsfile.c_str(),"a");
      int effect_fd = 0;
      effect_fd = fileno(effect_fh);
      if ( effect_fd == -1 ) 
	{ 
	  std::cerr << "ERROR: could not open " << params.effectsfile << '\n';
	  exit(10);
	}
      std::ostringstream indexstream;
      indexstream << params.replicate_no << ' ' << ftell(effect_fh) << ' '
		  << ftell(pheno_fh) << ' ' << ftell(haps_fh);
      fprintf(index_fh,"%s\n",indexstream.str().c_str());

      //write the haplotype data
      if ( ::write(hapfile_fd,popbuffer.str().c_str(),popbuffer.str().size()) == -1 )
	{
	  cerr << "Error writing to " << params.hapfile << '\n';
	  exit(errno);
	}
      //write the phenotype data
      if( ::write(pheno_fd, phenobuffer.str().c_str(), phenobuffer.str().size() ) == -1 )
	{
	  cerr << "Error writing to " << params.phenofile << '\n';
	  exit(errno);
	}
      //write the effects
      if( ::write( effect_fd, effectstream.str().c_str(), effectstream.str().size() ) == -1 )
	{
	  cerr << "Error writing to " << params.effectsfile << '\n';
	  exit(errno);
	}
      fflush( haps_fh );
      fclose( haps_fh );
      fflush( pheno_fh );
      fclose( pheno_fh );
      fflush( effect_fh );
      fclose( effect_fh );
    }

  
  // hapfile_flock.l_type = F_WRLCK;/*Write lock*/
  // hapfile_flock.l_whence = SEEK_SET;
  // hapfile_flock.l_start = 0;
  // hapfile_flock.l_len = 0;/*Lock whole file*/

  // phenofile_flock.l_type = F_WRLCK;/*Write lock*/
  // phenofile_flock.l_whence = SEEK_SET;
  // phenofile_flock.l_start = 0;
  // phenofile_flock.l_len = 0;/*Lock whole file*/

  // effects_flock.l_type = F_WRLCK;/*Write lock*/
  // effects_flock.l_whence = SEEK_SET;
  // effects_flock.l_start = 0;
  // effects_flock.l_len = 0;/*Lock whole file*/
  
  // FILE * haps_fh = fopen(params.hapfile.c_str(),"a");
  // int hapfile_fd = fileno(haps_fh);

  // if ( hapfile_fd == -1 ) 
  //   { 
  //     std::cerr << "ERROR: could not open " << params.hapfile << '\n';
  //     exit(10);
  //   }
  /*
  if (fcntl(index_fd, F_SETLKW,&hapfile_flock) == -1) 
    {
      std::cerr << "ERROR: could not obtain lock on " << params.hapfile << '\n';
      exit(10);
    }
  */
  // FILE * pheno_fh = fopen(params.phenofile.c_str(),"a");
  // int pheno_fd = fileno(pheno_fh);
  // if ( pheno_fd == -1 ) 
  //   { 
  //     std::cerr << "ERROR: could not open " << params.phenofile << '\n';
  //     exit(10);
  //   }
  /*
  if (fcntl(pheno_fd, F_SETLKW,&phenofile_flock) == -1) 
    {
      std::cerr << "ERROR: could not obtain lock on " << params.phenofile << '\n';
      exit(10);
    }
  */
  // FILE * effect_fh = NULL;
  // int effect_fd = 0;
  // effect_fh = fopen(params.effectsfile.c_str(),"a");
  // effect_fd = fileno(effect_fh);
  // if ( effect_fd == -1 ) 
  //   { 
  //     std::cerr << "ERROR: could not open " << params.effectsfile << '\n';
  //     exit(10);
  //   }
  /*
  if (fcntl(effect_fd, F_SETLKW,&effects_flock) == -1) 
    {
      std::cerr << "ERROR: could not obtain lock on " << params.effectsfile << '\n';
      exit(10);
    }
  */
  //Write the index data
  // std::ostringstream indexstream;
  // indexstream << params.replicate_no << ' ' << ftell(effect_fh) << ' '
  // 	      << ftell(pheno_fh) << ' ' << ftell(haps_fh);
  // fprintf(index_fh,"%s\n",indexstream.str().c_str());


  // //write the haplotype data
  // if ( ::write(hapfile_fd,popbuffer.str().c_str(),popbuffer.str().size()) == -1 )
  //   {
  //     cerr << "Error writing to " << params.hapfile << '\n';
  //     exit(errno);
  //   }
  //write the phenotype data
  // if( ::write(pheno_fd, phenobuffer.str().c_str(), phenobuffer.str().size() ) == -1 )
  //   {
  //     cerr << "Error writing to " << params.phenofile << '\n';
  //     exit(errno);
  //   }
  //write the effects
  // if( ::write( effect_fd, effectstream.str().c_str(), effectstream.str().size() ) == -1 )
  //   {
  //     cerr << "Error writing to " << params.effectsfile << '\n';
  //     exit(errno);
  //   }

  //release the locks
  index_flock.l_type = F_UNLCK;
  //hapfile_flock.l_type = F_UNLCK;
  //phenofile_flock.l_type = F_UNLCK;
  //effects_flock.l_type = F_UNLCK;

  /*
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
*/
  if (fcntl(index_fd, F_UNLCK,&index_flock) == -1) 
    {
      std::cerr << "ERROR: could not release lock on " << params.indexfile << '\n';
      exit(10);
    }
  fflush( index_fh );
  fclose(index_fh);
  
  exit(0);
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
    ("constant,C","Model constant effect size.  Otherwise, exponential distribution is used")
    ("multiplicative,m","Use multiplicative model of Risch and colleagues.  Default is Thornton, Foran & Long (2013) recessive model")
    ("popgen,F","Use popgen-like multiplicative model to calculate trait value")
    ("dominance,d",value<double>(&rv.dominance)->default_value(0.0),"Assign dominance for popgen-like model. Not used for any other model and will be ignored.")
    ("additive,a","Use additive model to calculate phenotype.  Default is Thornton, Foran & Long (2013) recessive model")
    ("indexfile,i",value<string>(&rv.indexfile)->default_value(string()),"Name of index file")
    ("popfile,p",value<string>(&rv.hapfile)->default_value(string()),"Name of output file for population")
    ("phenotypes,P",value<string>(&rv.phenofile)->default_value(string()),"Name of output file for phenotypes")
    ("effectsfile,E",value<string>(&rv.effectsfile)->default_value(string()),"Name of output file for effect sizes of causative mutations")
    ("seed,S",value<unsigned>(&rv.seed)->default_value(0),"Random number seed (unsigned integer)")
    ("optimum",value<double>(&rv.optimum)->default_value(0.),"At onset of exponential growth, change optimium value of phenotype.  Default is no change.")
    ("gzout","Output will be in compressed binary.  Default is plain binary.")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if(argc == 1 || vm.count("help"))
    {
      cerr << desc << '\n';
      exit(0);
    }

  if (vm.count("multiplicative"))
    {
      rv.model = MULTIPLICATIVE;
    }
  if (vm.count("additive"))
    {
      if( rv.model != GENE_RECESSIVE )
	{
	  cerr << "Error, it looks like multiple phenotype models have been chosen.  Please choose one!\n";
	  exit(10);
	}
      rv.model = GENE_ADDITIVE;
    }
  if( vm.count("popgen") )
    {
      if( rv.model != GENE_RECESSIVE )
	{
	  cerr << "Error, it looks like multiple phenotype models have been chosen.  Please choose one!\\
n";    
	  exit(10);
	}
      rv.model = POPGEN;
    }

  if(vm.count("gzout"))
    {
      rv.gzoutput=true;
    }

  if( vm.count("constant") )
    {
      rv.dist_effects = false;
    }
  if( rv.indexfile.empty() || rv.hapfile.empty() || rv.phenofile.empty() )
    {
      if (rv.indexfile.empty())
	{
	  cerr << "Error: index file name required.  Use -h to see options\n";
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
