/*
  The version of the Simons et al. eqn 
  used in the Hernandez paper http://biorxiv.org/content/biorxiv/early/2015/03/01/015917.full.pdf

  Note: they are not the same!!!!
 */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <zlib.h>
#include <fwdpp/diploid.hh>
#include <diseaseSims/mutation_with_age.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/program_options.hpp>
#include <gsl/gsl_randist.h>
#include <sys/stat.h>
#include <TFL_fitness_models.hpp>

using namespace std;
using namespace boost::accumulators;
using mean_acc = accumulator_set<double, stats<tag::mean> >;
using diploid_t = std::pair<glist::const_iterator,glist::const_iterator>;

struct vxv1params
{
  MODEL m;
  double h;
  std::string popfile;
  unsigned nreps;
};

vxv1params parse_argv( int argc, char ** argv );
  
/*
  Functions for TFL2013 and additive trait calculators.
  The simulation uses calculators that return pairs of doubles.
  The pair is the G,E components of trait value.
  Here, we just need the G, so we do this as a quick fix:

  (These should later be added right to gene_based_model.hpp!!!)
*/
struct TFL2013g
{
  using return_type = double;
  inline double operator()( const glist::const_iterator & g1,
			    const glist::const_iterator & g2) const
  {
    double e1 = std::accumulate( g1->smutations.begin(),
				 g1->smutations.end(),
				 0.,
				 [](const double & a,
				    const gtype::mutation_list_type_iterator & b) { return a + b->s; } );
    double e2 = std::accumulate( g2->smutations.begin(),
				 g2->smutations.end(),
				 0.,
				 [](const double & a,
				    const gtype::mutation_list_type_iterator & b) { return a + b->s; } );
    return (sqrt(e1*e2));
  }
};

struct additiveg
{
  using return_type = double;
  inline double operator()( const glist::const_iterator & g1,
			    const glist::const_iterator & g2) const
  {
    return std::accumulate( g1->smutations.begin(),
			    g1->smutations.end(),
			    0.,
			    [](const double & a,
			       const gtype::mutation_list_type_iterator & b) { return a + b->s; } )
      + std::accumulate( g2->smutations.begin(),
			 g2->smutations.end(),
			 0.,
			 [](const double & a,
			    const gtype::mutation_list_type_iterator & b) { return a + b->s; } );
  }
};

vector<double> getG( const dipvector & diploids,
		     const std::function<double(const glist::const_iterator &,
						const glist::const_iterator &)> dipG )
{
  vector<double> rv;
  for_each( diploids.begin(),diploids.end(),[&rv,&dipG](const diploid_t & __d ) { rv.push_back( dipG(__d.first,__d.second) ); } );
  return rv;
}

int main( int argc, char ** argv)
{
  const vxv1params pars = parse_argv(argc, argv);

  /*
    We now need to determine how to calculate the genetic 
    component of a diploid's trait value.

    The default will be the gene-based recessive model, aka TFL2013
  */
  std::function<double(const glist::const_iterator &,
		       const glist::const_iterator &)> dipG = std::bind(TFL2013g(),std::placeholders::_1,std::placeholders::_2);
  //handle user options
  switch( pars.m )
    {
    case MODEL::GENE_RECESSIVE:
      //default, do nothing
      break;
    case MODEL::GENE_ADDITIVE:
      dipG = std::bind(additiveg(),std::placeholders::_1,std::placeholders::_2);
      break;
    case MODEL::MULTIPLICATIVE:
      dipG = std::bind(multiplicative_phenotype(),std::placeholders::_1,std::placeholders::_2);
      break;
    case MODEL::POPGEN:
      dipG = std::bind(popgen_phenotype(),std::placeholders::_1,std::placeholders::_2,pars.h);
      break;
    case MODEL::EYREWALKER:
      cerr << "Eyre-Walker model not implemented yet\n";
      exit(EXIT_SUCCESS);
      break;
    }

  //Read each population in and process it
  gzFile gzin = gzopen( pars.popfile.c_str(),"rb" );
  for( unsigned i = 0 ; i < pars.nreps ; ++i )
    {
      mlist mutations;
      glist gametes;
      dipvector diploids;
      KTfwd::read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
      auto Gvals = getG(diploids,dipG);
    }
  gzclose(gzin);
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
    ("additive,a","Assume additive genetic model.  (Default is TFL2013)")
    ("multiplicative,m","Assume multiplicative genetic model.  (Default is TFL2013)")
    ("popgen,P","Assume multiplicative genetic model with dominance.  (Default is TFL2013)")
    ("dominance,d",value<double>(&rv.h)->default_value(0.0),"Dominance for multiplicative model with dominance")
    ("eyrewalker,e","Use Eyre-Walker 2010 model (not implemented yet!)")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);
  
  if(argc == 1 || vm.count("help"))
    {
      cerr << desc << '\n';
      exit(0);
    }

  if(!vm.count("popfile"))
    {
      cerr << "Error, no popfile specified.\n";
      exit(EXIT_SUCCESS);
    }
  else
    {
      struct stat buf;
      if (stat(rv.popfile.c_str(), &buf) != 0) {
	cerr << "Error: " << rv.popfile
	     << " does not exist.\n";
	exit(EXIT_FAILURE);
      }
    }
  if( !rv.nreps )
    {
      cerr << "Error, number of replicates == 0.\n";
      exit(EXIT_SUCCESS);
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
