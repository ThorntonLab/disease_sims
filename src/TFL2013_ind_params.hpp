#ifndef __TLF_IND_PARAMS_HPP__
#define __TLF_IND_PARAMS_HPP__

#include <string>
#include <iosfwd>
//fitness models
enum class MODEL { GENE_RECESSIVE = 0, GENE_ADDITIVE, MULTIPLICATIVE, POPGEN, EYREWALKER };

struct mut_model_params
{
  double mu_disease,mu_neutral,s,shape;
  unsigned N_ancestral;
  bool dist_effects;
  mut_model_params(void);
};

struct simparams
{
  mut_model_params mmp;
  unsigned N,N2,ngens_burnin,ngens_evolve,ngens_evolve_growth,replicate_no,seed;
  double littler,s,sd,sd_s,optimum,dominance;
  MODEL model;
  std::string indexfile, hapfile, phenofile, effectsfile ;
  bool verbose;
  simparams(void);
  std::ostream & print(std::ostream &) const;
};

std::ostream & operator<<(std::ostream &, const simparams &);
std::ostream & operator<<(std::ostream &, const MODEL &);

simparams parse_command_line(const int & argc,
			     char ** argv);

#endif
