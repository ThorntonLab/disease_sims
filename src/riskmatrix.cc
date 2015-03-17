/*
  Read in population, output a genotype matrix 
  for risk mutations.
 */

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
#include <cstdint>
#include <sstream>
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
using var_acc = accumulator_set<double, stats<tag::variance> >;
using diploid_t = std::pair<glist::const_iterator,glist::const_iterator>;
using Gfxn_t = std::function<double(const glist::const_iterator &,
				    const glist::const_iterator &)>;
struct riskmatrix_params
{
  MODEL m;
  double h;
  std::string popfile,outfile;
  unsigned nreps;
};

using vmcount_t = vector<pair<mlist::iterator,int8_t> >;

riskmatrix_params parse_argv( int argc, char ** argv );
Gfxn_t set_model( const riskmatrix_params & pars );
vector<double> getG( const dipvector & diploids,
		     const Gfxn_t & dipG );

vmcount_t get_mut_counts( const glist::const_iterator & g1,
			  const glist::const_iterator & g2 );

int main( int argc, char ** argv)
{
  const riskmatrix_params pars = parse_argv(argc, argv);

  /*
    We now need to determine how to calculate the genetic 
    component of a diploid's trait value.

    The default will be the gene-based recessive model, aka TFL2013
  */
  Gfxn_t dipG = set_model(pars);

  //Read each population in and process it
  gzFile gzin = gzopen( pars.popfile.c_str(),"rb" );
  for( unsigned rep = 0 ; rep < pars.nreps ; ++rep )
    {
      mlist mutations;
      glist gametes;
      dipvector diploids;
      KTfwd::read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
      unsigned RISKMUTIDX=0;
      vector<pair<mlist::iterator,unsigned> > risk_indexes;
      for( auto i = mutations.begin();i!=mutations.end();++i )
	{
	  if( ! i->neutral )
	    {
	      risk_indexes.push_back( make_pair(i,RISKMUTIDX++) );
	    }
	}
      auto Gvals = getG(diploids,dipG);
      //vector< vector<int8_t > > genos( diploids.size(), vector<int8_t>(RISKMUTIDX,0) );
      ostringstream buffer;
      for( unsigned ind = 0 ; ind < diploids.size() ; ++ind )
	{
	  vmcount_t vmc = get_mut_counts(diploids[ind].first,diploids[ind].second);
	  vector< int8_t > genos(RISKMUTIDX,0);
	  for( unsigned i = 0 ; i < vmc.size() ; ++i )
	    {
	      auto __itr = find_if( risk_indexes.begin(), risk_indexes.end(),[&vmc,&i](const pair<mlist::iterator,unsigned> & __p) {
		  return __p.first == vmc[i].first;
		});
	      genos[__itr->second] += vmc[i].second;
	      //genos[ind][distance(mutations.begin(),vmc[i].first)] += vmc[i].second;
	    }
	  buffer << Gvals[ind] << ' ';
	  for( unsigned i = 0 ; i < RISKMUTIDX ;++i )
	    {
	      buffer << int(genos[i]);
	      if(i < RISKMUTIDX-1) buffer << ' ';
	    }
	  buffer << '\n';
	}
      ostringstream ofname;
      ofname << pars.outfile << '.'
	     << rep << ".gz";
      gzFile gzout = gzopen(ofname.str().c_str(),"w");
      gzwrite(gzout,buffer.str().c_str(),buffer.str().size());
      gzclose(gzout);
    }
  gzclose(gzin);
}

riskmatrix_params parse_argv( int argc, char ** argv )
{
  using namespace boost::program_options;
  
  riskmatrix_params rv;
  rv.m = MODEL::GENE_RECESSIVE;
  rv.outfile = "riskmatrix_out";
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
    ("outfile,o",value<string>(&rv.outfile)->default_value("riskmatrix_out"),"Output file name base")
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

Gfxn_t set_model( const riskmatrix_params & pars )
{
  Gfxn_t dipG = std::bind(TFL2013g(),std::placeholders::_1,std::placeholders::_2);
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
  return dipG;
}

vector<double> getG( const dipvector & diploids,
		     const Gfxn_t & dipG )
{
  vector<double> rv;
  for_each( diploids.begin(),diploids.end(),[&rv,&dipG](const diploid_t & __d ) { rv.push_back( dipG(__d.first,__d.second) ); } );
  return rv;
}

vmcount_t get_mut_counts( const glist::const_iterator & g1,
			  const glist::const_iterator & g2 )
{
  vmcount_t rv;

  auto updater = [&rv](const mlist::iterator & __mut) {
    auto __itr =  find_if(rv.begin(),rv.end(),[&__mut](const pair<mlist::iterator,unsigned> & __p) {
	return __p.first == __mut;
      } );
    if(__itr == rv.end())
      rv.push_back(make_pair(__mut,1u));
    else
      __itr->second++;
  };
  for_each( g1->smutations.begin(), g1->smutations.end(),updater );
  for_each( g2->smutations.begin(), g2->smutations.end(),updater );
  return rv;
}

