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
using mean_acc = accumulator_set<double, stats<tag::mean> >;
using diploid_t = std::pair<glist::const_iterator,glist::const_iterator>;
using Gfxn_t = std::function<double(const glist::const_iterator &,
				    const glist::const_iterator &)>;
struct vxv1params
{
  MODEL m;
  double h;
  std::string popfile,outfile;
  unsigned nreps;
};

struct mphenos
{
  mean_acc AA,Aa,aa;
  mphenos() : AA(mean_acc()),Aa(mean_acc()),aa(mean_acc()) {}
};

struct running_means
{
  mean_acc z,d; //effect sizes, dominance
  unsigned nm; //number of mutations
  running_means() : z(mean_acc()),d(mean_acc()),nm(0u) {}
};

using mpheno_t = pair<mlist::iterator,mphenos>;
using mutphenovec_t = vector<mpheno_t>;
using vmcount_t = vector<pair<mlist::iterator,int8_t> >;

vxv1params parse_argv( int argc, char ** argv );
Gfxn_t set_model( const vxv1params & pars );
vector<double> getG( const dipvector & diploids,
		     const Gfxn_t & dipG );
mutphenovec_t init_mphenovec( mlist & mutations );
vmcount_t get_mut_counts( const glist::const_iterator & g1,
			  const glist::const_iterator & g2 );
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

int main( int argc, char ** argv)
{
  const vxv1params pars = parse_argv(argc, argv);

  /*
    We now need to determine how to calculate the genetic 
    component of a diploid's trait value.

    The default will be the gene-based recessive model, aka TFL2013
  */
  Gfxn_t dipG = set_model(pars);

  //Read each population in and process it
  gzFile gzin = gzopen( pars.popfile.c_str(),"rb" );
  map< double, running_means > data;
  for( unsigned rep = 0 ; rep < pars.nreps ; ++rep )
    {
      mlist mutations;
      glist gametes;
      dipvector diploids;
      KTfwd::read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
      auto Gvals = getG(diploids,dipG);
      mutphenovec_t mpv = init_mphenovec(mutations);
      for( auto & __m : mpv )
	{
	  for( unsigned ind = 0 ; ind < diploids.size() ; ++ind )
	    {
	      unsigned c = std::count(diploids[ind].first->smutations.begin(),
				      diploids[ind].first->smutations.end(),__m.first) +
		std::count(diploids[ind].second->smutations.begin(),
			   diploids[ind].second->smutations.end(),__m.first);
	      if( ! c )
		{
		  __m.second.aa( Gvals[ind] );
		}
	      else if (c == 1)
		{
		  __m.second.Aa( Gvals[ind] );
		}
	      else if (c == 2)
		{
		  __m.second.AA( Gvals[ind] );
		}
	      else
		{
		  cerr << "Error: impossible number of copies of a single mutation in a diploid ("
		       << c << "), line " << __LINE__ << " of " << __FILE__ << '\n';
		  exit(EXIT_FAILURE);
		}
	    }
	  double mAA = mean(__m.second.AA),mAa=mean(__m.second.Aa),maa=mean(__m.second.aa);
	  if( isfinite(mAA) && isfinite(mAa) && isfinite(maa) )
	    {
	      double dist = abs(mAA-maa);
	      double d2 = dist/2.;
	      //double d = 0.5*mAa/d2;
	      //double d = (mAa/mAA)*d2;
	      double pa = dist/2.,ma=-dist/2.;
	      double d = pa+(mAa - mAA);
	      double p = double(__m.first->n)/(2.*diploids.size());
	      double esize = pa + d*((1.-p)-p);

	      auto __ditr = data.find(p);
	      if(__ditr == data.end())
		{
		  data.insert(make_pair(p,running_means()));
		  __ditr=data.find(p);
		  __ditr->second.z(pow(esize,2.));
		  __ditr->second.d(d);
		  __ditr->second.nm++;
		}
	      else
		{
		  __ditr->second.z(pow(esize,2.));
		  __ditr->second.d(d);
		  __ditr->second.nm++;
		}
	      /*
	      cout << p << ' ' << __m.first->s << ' ' << mAA << ' ' << mAa << ' ' << maa << ' '
		   << pa/d2 << ' ' << d << ' ' << ma/d2 << ' ' << esize << ' '
		   << " | "
		   << d*mAA/d2 << '\n';
	      */
	    }
	}
    }
  gzclose(gzin);

  unsigned ttl_muts  = 0;
  for( auto ditr = data.begin() ; ditr != data.end() ; ++ditr ) ttl_muts += ditr->second.nm;

  ostringstream buffer;
  buffer << "p\tVp\td\n";
  double SUM=0.;
  for( auto ditr = data.begin() ; ditr != data.end() ; ++ditr )
    {
      double fx = double(ditr->second.nm)/double(ttl_muts);
      SUM += 0.5*mean(ditr->second.z)*fx*(ditr->first)*(1.-ditr->first);
      buffer << ditr->first << '\t' << SUM << '\t' << mean(ditr->second.d) << '\n';
    }
  gzin = gzopen(pars.outfile.c_str(),"w");
  gzwrite(gzin,buffer.str().c_str(),buffer.str().size());
  gzclose(gzin);
}

vxv1params parse_argv( int argc, char ** argv )
{
  using namespace boost::program_options;
  
  vxv1params rv;
  rv.m = MODEL::GENE_RECESSIVE;
  rv.outfile = "vxv1.out.gz";  
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

Gfxn_t set_model( const vxv1params & pars )
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

mutphenovec_t init_mphenovec( mlist & mutations )
{
  mutphenovec_t rv;
  for( auto i = mutations.begin();i!=mutations.end();++i )
    {
      rv.push_back(make_pair(i,mphenos()));
    }
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
