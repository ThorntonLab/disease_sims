#ifndef __TLF_IND_PARAMS_HPP__
#define __TLF_IND_PARAMS_HPP__

#include <string>
//fitness models
enum class MODEL { GENE_RECESSIVE = 0, GENE_ADDITIVE, MULTIPLICATIVE, POPGEN, EYREWALKER };

struct mut_model_params
{
  double mu_disease,mu_neutral,s,shape;
  unsigned * N_current;
  bool dist_effects;
  mut_model_params(void);
};

struct simparams
{
  mut_model_params mmp;
  unsigned N,N2,ngens_burnin,ngens_evolve,ngens_evolve_growth,replicate_no,seed;
  double littler,s,sd,sd_s,optimum,dominance;
  //double mu_disease,mu_neutral,littler,s,sd,sd_s,optimum,dominance;
  //bool dist_effects;
  MODEL model;
  std::string indexfile, hapfile, phenofile, effectsfile ;
  simparams(void);
};

simparams parse_command_line(const int & argc,
			     char ** argv);

#endif
