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

using namespace std;
using namespace boost::interprocess;
using namespace KTfwd;

//"tags" for the mutation model
struct mut_model_tag {};
struct fixed_tag : mut_model_tag {};
struct dist_tag : mut_model_tag {};
struct ew_tag : mut_model_tag {};  //Eyre-Walker 2010

//The mutation model -- old version
// struct mutation_model
// {
//   typedef TFLmtype result_type;
//   inline result_type operator()( gsl_rng * r, const unsigned int & ttl_generations,
//   				 const double & s, const double & ud, const double & un,
//   				 lookup_table_type * lookup,
//   				 const bool dist_effects = false) const
//   {
//     double pos = gsl_rng_uniform(r);
//     while(lookup->find(pos) != lookup->end())
//       {
// 	pos = gsl_rng_uniform(r);
//       }
//     lookup->insert(pos);
//     if( gsl_rng_uniform(r) <= ud/(ud+un) )
//       {
// 	if( ! dist_effects )
// 	  {
// 	    return TFLmtype(pos,s,1,ttl_generations,'A',false);
// 	  }
// 	else
// 	  {
// 	    return TFLmtype(pos,gsl_ran_exponential(r,s),1,ttl_generations,'A',false);
// 	  }
//       }
//     return TFLmtype(pos,0.,1,ttl_generations,'S',true);
//   }
// };

//New version, tag/dispatch model

inline double get_unique_pos(gsl_rng * r, lookup_table_type * lookup)
{
  double pos = gsl_rng_uniform(r);
  while(lookup->find(pos) != lookup->end())
    {
      pos = gsl_rng_uniform(r);
    }
  lookup->insert(pos);
  return pos;
}

template<typename tag_type> 
TFLmtype mut_model_details(gsl_rng * r, const unsigned int & ttl_generations,
			   const mut_model_params & mmp,
			   lookup_table_type * lookup)
{
  return TFLmtype();
}

template<>
inline TFLmtype mut_model_details<fixed_tag>(gsl_rng * r, const unsigned int & ttl_generations,
					     const mut_model_params & mmp,
					     lookup_table_type * lookup)
{
  double pos = get_unique_pos(r,lookup);
  if( gsl_rng_uniform(r) <= mmp.mu_disease/(mmp.mu_disease+mmp.mu_neutral) )
    {
      return TFLmtype(pos,mmp.s,1,ttl_generations,'A',false);      
    }
  return TFLmtype(pos,0.,1,ttl_generations,'S',true);
}

template<>
inline TFLmtype mut_model_details<dist_tag>(gsl_rng * r, const unsigned int & ttl_generations,
					    const mut_model_params & mmp,
					    lookup_table_type * lookup)
{
  double pos = get_unique_pos(r,lookup);
  if( gsl_rng_uniform(r) <= mmp.mu_disease/(mmp.mu_disease+mmp.mu_neutral) )
    {
      return TFLmtype(pos,gsl_ran_exponential(r,mmp.s),1,ttl_generations,'A',false);
    }
  return TFLmtype(pos,0.,1,ttl_generations,'S',true);
}

template<>
inline TFLmtype mut_model_details<ew_tag>(gsl_rng * r, const unsigned int & ttl_generations,
					  const mut_model_params & mmp,
					  lookup_table_type * lookup)
{
  double pos = get_unique_pos(r,lookup);
  if( gsl_rng_uniform(r) <= mmp.mu_disease/(mmp.mu_disease+mmp.mu_neutral) )
    {
      //4Ns ~ \Gamma with shape mmp.shape and mean mmp.s, so s = 4Ns/4N...
      double s = gsl_ran_gamma(r,mmp.shape,mmp.s/mmp.shape)/(4.*double(mmp.N_ancestral));
      return TFLmtype(pos,s,1,ttl_generations,'A',false);
    }
  return TFLmtype(pos,0.,1,ttl_generations,'S',true);
}

template<typename tag_type>
TFLmtype mut_model2(gsl_rng * r, const unsigned int & ttl_generations,
		    const mut_model_params & mmp,
		    lookup_table_type * lookup)
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

  lookup_table_type lookup;
  //the population begins with 1 gamete with no mutations
  glist gametes(1,gtype(2*params.N));
  mlist mutations;
  mvector fixations;      
  ftvector fixation_times;
  dipvector diploids(params.N,
		     std::make_pair(gametes.begin(),
				    gametes.begin()));
  unsigned generation;
  unsigned ttl_gen = 0;
  double wbar=1;

  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r);

  std::function<TFLmtype(mlist *)> MMODEL = std::bind(mut_model2<fixed_tag>,r,ttl_gen,params.mmp,&lookup);

  unsigned N_current = params.N;
  unsigned N_next = N_current;

  //Critical: mmp is getting bound to policies here, so we need to set the pointer to current N PRIOR to binding
  params.mmp.N_ancestral = params.N;
  if(params.mmp.dist_effects)
    MMODEL = std::bind(mut_model2<dist_tag>,r,ttl_gen,params.mmp,&lookup);
  if(params.model == MODEL::EYREWALKER)
    MMODEL = std::bind(mut_model2<ew_tag>,r,ttl_gen,params.mmp,&lookup);

  for( generation = 0; generation < params.ngens_burnin; ++generation,++ttl_gen )
    {
      //Evolution w/no deleterious mutations and no selection.
      wbar = sample_diploid(r,
			    &gametes,
			    &diploids,
			    &mutations,
			    params.N,
			    params.mmp.mu_neutral,
			    MMODEL,
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


  //Fitness model for phase w/selection.  The default is the recessive model of TFL 2013
  std::function<double(const glist::const_iterator &,
		       const glist::const_iterator &)> dipfit = std::bind(TFL2013_recessive(),std::placeholders::_1,std::placeholders::_2,params.sd,params.sd_s,0.,r);

  if( params.model == MODEL::GENE_ADDITIVE )
    {
      dipfit = std::bind(TFL2013_additive(),std::placeholders::_1,std::placeholders::_2,params.sd,params.sd_s,params.optimum,r);
    }
  else if( params.model == MODEL::MULTIPLICATIVE )
    {
      dipfit = std::bind(multiplicative_disease_effect_to_fitness(),std::placeholders::_1,std::placeholders::_2,
			 params.sd,params.sd_s,params.optimum,r);
    }
  else if ( params.model == MODEL::POPGEN )
    {
      dipfit = std::bind(popgen_disease_effect_to_fitness(),std::placeholders::_1,std::placeholders::_2,params.dominance,
			 params.sd,params.sd_s,params.optimum,r);
    }
  else if ( params.model == MODEL::EYREWALKER )
    {
      dipfit = std::bind(EWfitness(),std::placeholders::_1,std::placeholders::_2);
    }

  for( generation = 0; generation < params.ngens_evolve; ++generation,++ttl_gen )
    {
      //Evolve under the disease model from TFL2013
      wbar = sample_diploid(r,
			    &gametes,
			    &diploids,
			    &mutations,
			    params.N,
			    params.mmp.mu_disease+params.mmp.mu_neutral,
			    MMODEL,
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
      N_next = round( params.N*pow(G,generation+1) );
      wbar = sample_diploid(r,
			    &gametes,
			    &diploids,
			    &mutations,
			    N_current,
			    N_next,
			    params.mmp.mu_disease+params.mmp.mu_neutral,
			    MMODEL,
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
  if (params.model != MODEL::EYREWALKER) 
    {
      phenobuffer.write( reinterpret_cast<char *>(&nphenos), sizeof(unsigned) );
      for( unsigned i = 0 ; i < diploids.size() ; ++i )
	{
	  if (params.model == MODEL::GENE_RECESSIVE) 
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
	  if (params.model == MODEL::GENE_ADDITIVE) 
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
	  else if (params.model == MODEL::MULTIPLICATIVE) 
	    {
	      double x = multiplicative_phenotype()(diploids[i].first,
						    diploids[i].second);
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	      x = (params.sd > 0.) ? gsl_ran_gaussian(r,params.sd) : 0.;
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	    }
	  else if (params.model == MODEL::POPGEN)
	    {
	      double x = popgen_phenotype()(diploids[i].first,
					    diploids[i].second,params.dominance);
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	      x = (params.sd > 0.) ? gsl_ran_gaussian(r,params.sd) : 0.;
	      phenobuffer.write( reinterpret_cast< char * >(&x), sizeof(double) );
	    }
	}
    }

  //Write out effects information for causative sites
  ostringstream effectstream;
  unsigned nmuts = mutations.size();

  //effectstream.write( reinterpret_cast<char *>(&ncausative),sizeof(unsigned) );
  effectstream.write( reinterpret_cast<char *>(&nmuts),sizeof(unsigned) );
  for( typename mlist::const_iterator i = mutations.begin() ; i != mutations.end() ; ++i )
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
  //lock index file ASAP
  ofstream indexstream(params.indexfile.c_str(),ios::out|ios::app);
  file_lock flock(params.indexfile.c_str());
  scoped_lock<file_lock> slock(flock);

  gzFile gzout = gzopen(params.hapfile.c_str(),"a");
  int hapswritten = gzwrite(gzout,popbuffer.str().c_str(),popbuffer.str().size());
  if ( ! hapswritten && !popbuffer.str().empty() )
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
