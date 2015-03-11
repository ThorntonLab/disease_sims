/*
  The version of the Simons et al. eqn 
  used in the Hernandez paper http://biorxiv.org/content/biorxiv/early/2015/03/01/015917.full.pdf

  Note: they are not the same!!!!
 */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <map>
#include <zlib.h>
#include <fwdpp/diploid.hh>
#include <diseaseSims/mutation_with_age.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/program_options.hpp>
#include <gsl/gsl_randist.h>

#include <TFL_fitness_models.hpp>

using namespace std;
using namespace boost::accumulators;

struct vxv1params
{
  MODEL m;
  double h;
  std::string popfile;
  unsigned nreps;
};

vxv1params parse_argv( int argc, char ** argv );
  
using mean_acc = accumulator_set<double, stats<tag::mean> >;

int main( int argc, char ** argv)
{
  const vxv1params pars = parse_argv(argc, argv);
}

vxv1params parse_argv( int argc, char ** argv )
{
  using namespace boost::program_options;
  
  vxv1params rv;
  rv.m = MODEL::GENE_RECESSIVE;
  
  options_description desc("Calculate stuff");
  desc.add_options()
    ("help,h","Print usage information to screen")
    ("popfile,p",value<string>(&rv.popfile),"Name of file containing simulation output")
    ("nreps,n",value<unsigned>(&rv.nreps)->default_value(0u),"Number of replicates stored in popfile")
    ("additive","Assume additive genetic model.  (Default is TFL2013)")
    ("multiplicative","Assume multiplicative genetic model.  (Default is TFL2013)")
    ("popgen","Assume multiplicative genetic model with dominance.  (Default is TFL2013)")
    ("dominance",value<double>(&rv.h)->default_value(0.0),"Dominance for multiplicative model with dominance")
    ("eyrewalker","Use Eyre-Walker 2010 model (not implemented yet!)")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);
  
  if(argc == 1 || vm.count("help"))
    {
      cerr << desc << '\n';
      exit(0);
    }

  if(vm.count("additive"))
    {
      rv.m = MODEL::GENE_ADDITIVE;
    }
  else if(vm.count("multiplicative"))
    {
      rv.m = MODEL::MULTIPLICATIVE;
    }
  else if(vm.count("popgen"))
    {
      rv.m = MODEL::POPGEN;
    }

  if(vm.count("eyrewalker"))
    {
      cerr << "Eyre-Walker 2010 not implemented yet!\n";
      exit(EXIT_SUCCESS);
    }
  return rv;
}
