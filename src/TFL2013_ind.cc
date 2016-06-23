/*
  Simulation from Thornton, Foran, and Long (2013) PLoS Genetics.

  Re-implemented using individual-based sampling routines.

  Rewritten to use the fwdpp library v >= 0.2.4
*/
#include <fwdpp/diploid.hh>
#include <utility>
#include <iostream>
#include <fstream>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>

#include <zlib.h>

#include <TFL2013_ind_params.hpp>
#include <diseaseSims/mutation_with_age.hpp>
#include <TFL_fitness_models.hpp>

#include <fwdpp/sugar/serialization.hpp>
#include <fwdpp/sugar/infsites.hpp>
using namespace std;
using namespace boost::interprocess;
using namespace KTfwd;

//"tags" for the mutation model
struct mut_model_tag {};
struct fixed_tag : mut_model_tag {};
struct dist_tag : mut_model_tag {};
struct ew_tag : mut_model_tag {};  //Eyre-Walker 2010

//New version, tag/dispatch model

inline double get_unique_pos(gsl_rng * r, poptype::lookup_table_t * lookup)
{
  double pos = gsl_rng_uniform(r);
  while(lookup->find(pos) != lookup->end())
    {
      pos = gsl_rng_uniform(r);
    }
  lookup->insert(pos);
  return pos;
}

/*
  This is a nonsensical version of the function. Essentially a place-holder
  for the explicit speclizations that come later.
*/
template<typename tag_type>
TFLmtype mut_model_details(gsl_rng * r, const unsigned int & ttl_generations,
			   const mut_model_params & mmp,
			   poptype::lookup_table_t * lookup)
{
  //Return a nonsensical value that will definitely cause assertions to fail, etc. :)
  return TFLmtype(numeric_limits<double>::quiet_NaN(),
		  numeric_limits<double>::quiet_NaN(),
		  numeric_limits<double>::quiet_NaN(),
		  numeric_limits<unsigned>::max());
}

template<>
inline TFLmtype mut_model_details<fixed_tag>(gsl_rng * r, const unsigned int & ttl_generations,
					     const mut_model_params & mmp,
					     poptype::lookup_table_t * lookup)
{
  double pos = get_unique_pos(r,lookup);
  if( gsl_rng_uniform(r) <= mmp.mu_disease/(mmp.mu_disease+mmp.mu_neutral) )
    {
      return TFLmtype(pos,mmp.s,1.0,ttl_generations);//,'A',false);      
    }
  return TFLmtype(pos,0.,1.0,ttl_generations);//,'S',true);
}

template<>
inline TFLmtype mut_model_details<dist_tag>(gsl_rng * r, const unsigned int & ttl_generations,
					    const mut_model_params & mmp,
					    poptype::lookup_table_t * lookup)
{
  double pos = get_unique_pos(r,lookup);
  if( gsl_rng_uniform(r) <= mmp.mu_disease/(mmp.mu_disease+mmp.mu_neutral) )
    {
      //return TFLmtype(pos,gsl_ran_exponential(r,mmp.s),ttl_generations,'A',false);
      return TFLmtype(pos,gsl_ran_exponential(r,mmp.s),1.0,ttl_generations);//,'A',false);
    }
  return TFLmtype(pos,0.,1.0,ttl_generations);//,'S',true);
}

template<>
inline TFLmtype mut_model_details<ew_tag>(gsl_rng * r, const unsigned int & ttl_generations,
					  const mut_model_params & mmp,
					  poptype::lookup_table_t * lookup)
{
  double pos = get_unique_pos(r,lookup);
  if( gsl_rng_uniform(r) <= mmp.mu_disease/(mmp.mu_disease+mmp.mu_neutral) )
    {
      //4Ns ~ \Gamma with shape mmp.shape and mean mmp.s, so s = 4Ns/4N...
      double s = gsl_ran_gamma(r,mmp.shape,mmp.s/mmp.shape)/(4.*double(mmp.N_ancestral));
      return TFLmtype(pos,s,1.0,ttl_generations);//,'A',false);
      //return TFLmtype(pos,s,ttl_generations,'A',false);
    }
  return TFLmtype(pos,0.,1.0,ttl_generations);
  //return TFLmtype(pos,0.,ttl_generations,'S',true);
}

template<typename tag_type>
TFLmtype mut_model2(gsl_rng * r, const unsigned int & ttl_generations,
		    const mut_model_params & mmp,
		    poptype::lookup_table_t * lookup)
{
  static_assert(std::is_base_of<mut_model_tag,tag_type>::value,"tag_type must be derived from mut_model_tag");
  return mut_model_details<tag_type>(r,ttl_generations,mmp,lookup);
}

int main(int argc, char ** argv)
{
  simparams params = parse_command_line(argc,argv);
  if( params.verbose )
    {
      cerr << params << '\n';
      exit(EXIT_SUCCESS);
    }

  //Determine growth rate under exponential growth model.
  double G = 0.;
  if( params.ngens_evolve_growth > 0 ) 
    { 
      G = exp( (log(double(params.N2)) - log(double(params.N)))/double(params.ngens_evolve_growth)); 
    }

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,params.seed);

  poptype pop(2*params.N);
  unsigned generation;
  unsigned ttl_gen = 0;
  double wbar=1;

  //std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r);

  //std::function<TFLmtype(poptype::mcont_t &)> MMODEL = std::bind(mut_model2<fixed_tag>,r,ttl_gen,params.mmp,&pop.mut_lookup);
  
  //mutation model with constant effect size...
  std::function<std::size_t(std::queue<std::size_t> &,poptype::mcont_t &)> MMODEL = std::bind(KTfwd::infsites(),
											      std::placeholders::_1,std::placeholders::_2,r,
											      std::ref(pop.mut_lookup),&generation,
											      params.mmp.mu_neutral,
											      params.mmp.mu_disease,
											      [&r](){return gsl_rng_uniform(r);},
											      [&params](){return params.mmp.s;},
											      [](){return 0.;});
  unsigned N_current = params.N;
  unsigned N_next = N_current;

  //Critical: mmp is getting bound to policies here, so we need to set the pointer to current N PRIOR to binding
  params.mmp.N_ancestral = params.N;
  if(params.mmp.dist_effects)
    MMODEL = std::bind(KTfwd::infsites(),
		       std::placeholders::_1,std::placeholders::_2,r,
		       std::ref(pop.mut_lookup),&generation,
		       params.mmp.mu_neutral,
		       params.mmp.mu_disease,
		       [&r](){return gsl_rng_uniform(r);},
		       [&r,&params](){return gsl_ran_exponential(r,params.mmp.s);},
		       [](){return 0.;});
  if(params.model == MODEL::EYREWALKER)
    MMODEL = std::bind(KTfwd::infsites(),
		       std::placeholders::_1,std::placeholders::_2,r,
		       std::ref(pop.mut_lookup),&generation,
		       params.mmp.mu_neutral,
		       params.mmp.mu_disease,
		       [&r](){return gsl_rng_uniform(r);},
		       [&r,&params](){return gsl_ran_gamma(r,params.mmp.shape,params.mmp.s/params.mmp.shape)/(4.*double(params.mmp.N_ancestral));},
		       [](){return 0.;});
    //double s = 
    //MMODEL = std::bind(mut_model2<ew_tag>,r,ttl_gen,params.mmp,&pop.mut_lookup);
  for( generation = 0; generation < params.ngens_burnin; ++generation,++ttl_gen )
    {
      //Evolution w/no deleterious mutations and no selection.
      wbar = sample_diploid(r,
			    pop.gametes,
			    pop.diploids,
			    pop.mutations,
			    pop.mcounts,
			    params.N,
			    params.mmp.mu_neutral,
			    MMODEL,
			    std::bind(KTfwd::poisson_xover(),r,params.littler,0.,1.,
				      std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
			    std::bind(KTfwd::no_selection(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
			    pop.neutral,pop.selected);
      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,2*params.N);
    }


  //Fitness model for phase w/selection.  The default is the recessive model of TFL 2013
  std::function<double(const poptype::gamete_t &,
		       const poptype::gamete_t &,
		       const poptype::mcont_t &)> dipfit = std::bind(TFL2013_recessive(),
								     std::placeholders::_1,
								     std::placeholders::_2,
								     std::placeholders::_3,
								     params.sd,params.sd_s,0.,r);

  if( params.model == MODEL::GENE_ADDITIVE )
    {
      dipfit = std::bind(TFL2013_additive(),std::placeholders::_1,
			 std::placeholders::_2,
			 std::placeholders::_3,
			 params.sd,params.sd_s,params.optimum,r);
    }
  else if( params.model == MODEL::MULTIPLICATIVE )
    {
      dipfit = std::bind(multiplicative_disease_effect_to_fitness(),
			 std::placeholders::_1,
			 std::placeholders::_2,
			 std::placeholders::_3,
			 params.sd,params.sd_s,params.optimum,r);
    }
  else if ( params.model == MODEL::POPGEN )
    {
      dipfit = std::bind(popgen_disease_effect_to_fitness(),
			 std::placeholders::_1,
			 std::placeholders::_2,
			 std::placeholders::_3,
			 params.dominance,
			 params.sd,params.sd_s,params.optimum,r);
    }
  else if ( params.model == MODEL::EYREWALKER )
    {
      dipfit = std::bind(EWfitness(),std::placeholders::_1,
			 std::placeholders::_2,
			 std::placeholders::_3);
    }

  for( generation = 0; generation < params.ngens_evolve; ++generation,++ttl_gen )
    {
      //Evolve under the disease model from TFL2013
      wbar = sample_diploid(r,
			    pop.gametes,
			    pop.diploids,
			    pop.mutations,
			    pop.mcounts,
			    params.N,
			    params.mmp.mu_disease+params.mmp.mu_neutral,
			    MMODEL,
			    std::bind(KTfwd::poisson_xover(),r,params.littler,0.,1.,
				      std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
			    dipfit,
			    pop.neutral,pop.selected);
      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,2*params.N);
    }
  //Exp. growth phase w/disease model
  for( generation = 0 ; generation < params.ngens_evolve_growth ; ++generation,++ttl_gen )
    {
      N_next = round( params.N*pow(G,generation+1) );
      wbar = sample_diploid(r,
			    pop.gametes,
			    pop.diploids,
			    pop.mutations,
			    pop.mcounts,
			    N_current,
			    N_next,
			    params.mmp.mu_disease+params.mmp.mu_neutral,
			    MMODEL,
			    std::bind(KTfwd::poisson_xover(),r,params.littler,0.,1.,
				      std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
			    dipfit,
			    pop.neutral,pop.selected);
      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,2*params.N);
      //update N
      N_current = N_next;
    }

  //Write out the population
  KTfwd::serialize serializer;
  serializer(pop,std::bind(KTfwd::mutation_writer(),std::placeholders::_1,std::placeholders::_2));
  //ostringstream popbuffer;
  //write_binary_pop(&gametes,&mutations,&diploids,std::bind(mwriter(),std::placeholders::_1,std::placeholders::_2),popbuffer);

  //Write out the phenotypes
  ostringstream phenobuffer;
  unsigned nphenos = pop.diploids.size();
  if (params.model != MODEL::EYREWALKER) 
    {
      phenobuffer.write( reinterpret_cast<char *>(&nphenos), sizeof(unsigned) );
      for( unsigned i = 0 ; i < pop.diploids.size() ; ++i )
	{
	  if (params.model == MODEL::GENE_RECESSIVE) 
	    {
	      pair<double,double> pheno = TFL2013_recessive_disease_effect()(pop.gametes[pop.diploids[i].first],
									     pop.gametes[pop.diploids[i].second],
									     pop.mutations,
									     params.sd,
									     r);
	      double x = pheno.first;
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	      x = pheno.second;
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	    }
	  if (params.model == MODEL::GENE_ADDITIVE) 
	    {
	      pair<double,double> pheno = TFL2013_additive_disease_effect()(pop.gametes[pop.diploids[i].first],
									     pop.gametes[pop.diploids[i].second],
									     pop.mutations,
									    params.sd,
									    r);
	      double x = pheno.first;
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	      x = pheno.second;
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	    }
	  else if (params.model == MODEL::MULTIPLICATIVE) 
	    {
	      double x = multiplicative_phenotype()(pop.gametes[pop.diploids[i].first],
						    pop.gametes[pop.diploids[i].second],
						    pop.mutations);
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	      x = (params.sd > 0.) ? gsl_ran_gaussian(r,params.sd) : 0.;
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	    }
	  else if (params.model == MODEL::POPGEN)
	    {
	      double x = popgen_phenotype()(pop.gametes[pop.diploids[i].first],
									     pop.gametes[pop.diploids[i].second],
									     pop.mutations,params.dominance);
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	      x = (params.sd > 0.) ? gsl_ran_gaussian(r,params.sd) : 0.;
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	    }
	}
    }

  //Write out effects information for causative sites
  ostringstream effectstream;
  unsigned nmuts = pop.mutations.size();

  //effectstream.write( reinterpret_cast<char *>(&ncausative),sizeof(unsigned) );
  effectstream.write( reinterpret_cast<char *>(&nmuts),sizeof(unsigned) );
  for( auto i = pop.mutations.begin() ; i != pop.mutations.end() ; ++i )
    {
      if(! params.effectsfile.empty() )
	{
	  //position of causative mutation, effect size, count, age
	  double pos = i->pos;
	  double s = i->s;
	  double count = double(pop.mcounts[i-pop.mutations.begin()]);
	  double age = double(ttl_gen - i->g + 1);
	  effectstream.write( reinterpret_cast<char *>(&pos), sizeof(double) );
	  effectstream.write( reinterpret_cast<char *>(&s), sizeof(double) );
	  effectstream.write( reinterpret_cast<char *>(&count), sizeof(double) );
	  effectstream.write( reinterpret_cast<char *>(&age), sizeof(double) );
	}
    }
  //lock index file ASAP
  ofstream indexstream(params.indexfile.c_str(),ios::out|ios::app);
  file_lock flock(params.indexfile.c_str());
  scoped_lock<file_lock> slock(flock);

  gzFile gzout = gzopen(params.hapfile.c_str(),"a");
  int hapswritten = gzwrite(gzout,serializer.buffer.str().c_str(),serializer.buffer.str().size());
  if ( ! hapswritten && !serializer.buffer.str().empty() )
    {
      cerr << "Error writing population to " 
	   << params.hapfile << '\n';
      exit(EXIT_FAILURE);
    }
  gzclose(gzout);

  int phenowritten = -1;
  if( params.model != MODEL::EYREWALKER )
    {
      gzout = gzopen(params.phenofile.c_str(),"a");
      phenowritten = gzwrite(gzout,phenobuffer.str().c_str(),phenobuffer.str().size());
      if ( ! phenowritten && !phenobuffer.str().empty() )
	{
	  cerr << "Error writing population to " 
	       << params.phenofile << '\n';
	  exit(EXIT_FAILURE);
	}
      gzclose(gzout);
    }


  gzout = gzopen(params.effectsfile.c_str(),"a");
  int effectwritten = gzwrite(gzout,effectstream.str().c_str(),effectstream.str().size());
  if ( ! effectwritten && !effectstream.str().empty() )
    {
      cerr << "Error writing population to " 
	   << params.effectsfile << '\n';
      exit(EXIT_FAILURE);
    }
  gzclose(gzout);
    

  //Now we can write to the index file
  indexstream << params.replicate_no << ' ' << effectwritten << ' '
	      << phenowritten << ' ' << hapswritten << '\n';
  //release the locks
  indexstream.flush();
  flock.unlock();
  indexstream.close();
  exit(0);
}
