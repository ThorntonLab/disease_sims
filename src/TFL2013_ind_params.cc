#include <TFL2013_ind_params.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <limits>
using namespace std;
using namespace boost::program_options;

mut_model_params::mut_model_params(void) : mu_disease(0.000125),
					   mu_neutral(0.00125),
					   s(0.1),
					   shape(std::numeric_limits<double>::min()),
					   N_ancestral(0u),
					   dist_effects(true)
{
}

simparams::simparams(void) : mmp(mut_model_params()),
			     N(20000),N2(20000),
			     ngens_burnin(0),
			     ngens_evolve(160000),
			     ngens_evolve_growth(0),
			     replicate_no(0),
			     seed(0),
			     littler(0.00125),
			     sd(0.075),
			     sd_s(1),
			     optimum(0.),
			     dominance(0.),
			     model( MODEL::GENE_RECESSIVE ),
			     indexfile(string()),
			     hapfile(string()),
			     phenofile(string()),
			     effectsfile(string())
{
}

void param_error(const char * param,
		 const char * condition,
		 const int & val)
{
  cerr << "Error: " << param << ' '
       << condition << '\n';
  exit(val);
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
    ("neutral,n",value<double>(&rv.mmp.mu_neutral)->default_value(0.00125),"Neutral mutation rate (per region per generation")
    ("causative,c",value<double>(&rv.mmp.mu_disease)->default_value(0.000125),"Mutation rate to causative mutations(per region per generation")
    ("recrate,r",value<double>(&rv.littler)->default_value(0.00125),"Recombination rate (per diploid per region per generation")
    ("esize,e",value<double>(&rv.mmp.s)->default_value(0.1),"Effect size of causative mutation.  Mean of exponential dist by default.  Constant effect size if dist 0 or -d 0 is used")
    ("noise",value<double>(&rv.sd)->default_value(0.075),"Std. deviation in Gaussian noise to be added to phenotype")
    ("sigma",value<double>(&rv.sd_s)->default_value(1.0),"Std. deviation in Gaussian fitness function")
    ("constant,C","Model constant effect size.  Otherwise, exponential distribution is used")
    ("multiplicative,m","Use multiplicative model of Risch and colleagues.  Default is Thornton, Foran & Long (2013) recessive model")
    ("popgen,F","Use popgen-like multiplicative model to calculate trait value")
    ("eyre-walker","Use the model froim Eyre-Walker 2010, www.pnas.org/cgi/doi/10.1073/pnas.0906182107")
    ("ewshape",value<double>(&rv.mmp.shape),"Shape parameter for Gamma density of S = 4Ns.  The value for --esize/-e will be taken as the mean value for S.")
    ("dominance,d",value<double>(&rv.dominance)->default_value(0.0),"Assign dominance for popgen-like model. Not used for any other model and will be ignored.")
    ("additive,a","Use additive model to calculate phenotype.  Default is Thornton, Foran & Long (2013) recessive model")
    ("indexfile,i",value<string>(&rv.indexfile)->default_value(string()),"Name of index file")
    ("popfile,p",value<string>(&rv.hapfile)->default_value(string()),"Name of output file for population")
    ("phenotypes,P",value<string>(&rv.phenofile)->default_value(string()),"Name of output file for phenotypes.  Not used for Eyre-Walker model.")
    ("effectsfile,E",value<string>(&rv.effectsfile)->default_value(string()),"Name of output file for effect sizes of causative mutations")
    ("seed,S",value<unsigned>(&rv.seed)->default_value(0),"Random number seed (unsigned integer)")
    ("optimum",value<double>(&rv.optimum)->default_value(0.),"At onset of exponential growth, change optimium value of phenotype.  Default is no change.  Not allowed for the Eyre-Walker model")
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
      rv.model = MODEL::MULTIPLICATIVE;
    }
  if (vm.count("additive"))
    {
      if( rv.model != MODEL::GENE_RECESSIVE )
	{
	  cerr << "Error, it looks like multiple phenotype models have been chosen.  Please choose one!\n";
	  exit(EXIT_FAILURE);
	}
      rv.model = MODEL::GENE_ADDITIVE;
    }
  if( vm.count("popgen") )
    {
      if( rv.model != MODEL::GENE_RECESSIVE )
	{
	  cerr << "Error, it looks like multiple phenotype models have been chosen.  Please choose one!\n";    
	  exit(EXIT_FAILURE);
	}
      rv.model = MODEL::POPGEN;
    }
  if( vm.count("eyre-walker") ) 
    {
      if( rv.model != MODEL::GENE_RECESSIVE )
	{
	  cerr << "Error, it looks like multiple phenotype models have been chosen.  Please choose one!\n";    
	  exit(EXIT_FAILURE);
	}
      rv.model = MODEL::EYREWALKER;
      if (! vm.count("ewshape") )
	{
	  cerr << "Error: shape parameter required for Eyre-Walker model.\n";
	  exit(EXIT_FAILURE);
	}
      else if (	rv.mmp.shape <= 0. )
	{
	  param_error("shape parameter for Gamma distribution (--ewshape) ","is <= 0",EXIT_FAILURE);
	}
      if(rv.optimum != 0.)
	{
	  cerr << "Error: optimum shift currently not allowed for the Eyre-Walker model.\n";
	  exit(EXIT_FAILURE);
	}
    }

  if( vm.count("constant") )
    {
      rv.mmp.dist_effects = false;
    }
  if( rv.indexfile.empty() || rv.hapfile.empty() || (rv.phenofile.empty()&&rv.model!=MODEL::EYREWALKER) )
    {
      if (rv.indexfile.empty())
	{
	  cerr << "Error: index file name required.  Use -h to see options\n";
	}

      if( rv.hapfile.empty() )
	{
	  cerr << "Error: population output file name required.  Use -h to see options\n";
	}
	
      if( rv.phenofile.empty() && rv.model != MODEL::EYREWALKER )
	{
	  cerr << "Error: phenotypes output file name required.  Use -h to see options\n";
	}
      exit(EXIT_FAILURE);
    }

  //Sanity checks on parameter values
  if( rv.mmp.mu_neutral < 0. )
    {
      param_error("neutral mutation rate (-n/--neutral)","is < 0",EXIT_FAILURE);
    }
  if( rv.mmp.mu_disease < 0. )
    {
      param_error("causative mutation rate (-c/--causative)","is < 0",EXIT_FAILURE);
    }
  if( rv.mmp.s < 0. )
    {
      param_error("effect size (-e/--esize)","is < 0",EXIT_FAILURE);
    }
  if( rv.sd < 0. )
    {
      param_error("Gaussian noise sigma (--noise)","is < 0",EXIT_FAILURE);
    }
  if( rv.sd_s < 0. )
    {
      param_error("Gaussian fitness function sigma (--sigma)","is < 0",EXIT_FAILURE);
    }
  return rv;
}
