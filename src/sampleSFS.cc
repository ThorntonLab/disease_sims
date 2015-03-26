#include <algorithm>
#include <functional>
#include <random>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <zlib.h>
#include <fwdpp/IO.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/program_options.hpp>
#include <diseaseSims/mutation_with_age.hpp>

using namespace std;
using namespace KTfwd;
using namespace boost::accumulators;

using mean_acc = accumulator_set<double, stats<tag::mean> >;
struct params 
{
  string infile;
  string outfile;
  unsigned nsam,seed;
};

params parse_argv(int argc, char ** argv);

int main(int argc, char ** argv)
{
  params pars = parse_argv(argc,argv);

  gzFile gzin = gzopen( pars.infile.c_str(),"rb" );
  if( gzin == NULL )
    {
      cerr << "Error, " << pars.infile 
	   << " could not be opened for reading in binary mode. "
	   << " Line " << __LINE__ << " of " << __FILE__ << '\n';
      exit(EXIT_FAILURE);
    }

  mt19937 generator;
  generator.seed(pars.seed);
  vector< mean_acc > meanRiskSFS(pars.nsam),meanNeutSFS(pars.nsam);

  auto updater = []( const mlist::iterator & __m, vector< pair<const mlist::iterator,unsigned> > & __counts  ) {
    auto __i = find_if(__counts.begin(),__counts.end(),[&__m](const pair<mlist::iterator,unsigned> & __p ) {
	return __m == __p.first;
      });
    if(__i == __counts.end())
      {
	__counts.push_back(make_pair(__m,1u));
      }
    else
      {
	__i->second++;
      }
  };

  do 
    {
      mlist mutations;
      glist gametes;
      dipvector diploids;
      KTfwd::read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
      vector<unsigned> rd(diploids.size());
      unsigned dummy = 0;
      for_each(rd.begin(),rd.end(),[&dummy](unsigned & __u){__u=dummy++;});
      shuffle(rd.begin(),rd.end(),generator);
      vector< pair<const mlist::iterator,unsigned> >  sampleRiskCounts,sampleNeutralCounts;

      for( unsigned i = 0 ; i < pars.nsam/2 ; ++i )
	{
	  assert( rd[i] < diploids.size() );
	  for_each( diploids[ rd[i] ].first->mutations.begin(),diploids[ rd[i] ].first->mutations.end(),std::bind(updater,std::placeholders::_1,ref(sampleNeutralCounts)) );
	  for_each( diploids[ rd[i] ].second->mutations.begin(),diploids[ rd[i] ].second->mutations.end(),std::bind(updater,std::placeholders::_1,ref(sampleNeutralCounts)) );
	  for_each( diploids[ rd[i] ].first->smutations.begin(),diploids[ rd[i] ].first->smutations.end(),std::bind(updater,std::placeholders::_1,ref(sampleRiskCounts)) );
	  for_each( diploids[ rd[i] ].second->smutations.begin(),diploids[ rd[i] ].second->smutations.end(),std::bind(updater,std::placeholders::_1,ref(sampleRiskCounts)) );
	}

      std::vector<unsigned> nSFS(pars.nsam,0u),rSFS(pars.nsam,0u);
      for( unsigned i = 0 ; i < sampleNeutralCounts.size() ; ++i )
	{
	  assert( sampleNeutralCounts[i].second > 0 );
	  assert( sampleNeutralCounts[i].second <= pars.nsam );
	  nSFS[ sampleNeutralCounts[i].second - 1 ]++;
	}
      for( unsigned i = 0 ; i < sampleRiskCounts.size() ; ++i )
	{
	  assert( sampleRiskCounts[i].second > 0 );
	  assert( sampleRiskCounts[i].second <= pars.nsam );
	  rSFS[ sampleRiskCounts[i].second - 1 ]++;
	}

      for( unsigned i = 0 ; i < nSFS.size() ; ++i ) meanNeutSFS[i](nSFS[i]);
      for( unsigned i = 0 ; i < rSFS.size() ; ++i ) meanRiskSFS[i](rSFS[i]);

      //Check if at EOF
      int c = gzgetc(gzin);
      if(gzeof(gzin)) break;
      else gzungetc(c,gzin);
    }
  while(! gzeof(gzin) );
  gzclose(gzin);

  ostringstream buffer;
  for ( unsigned i = 0 ; i < pars.nsam ; ++i )
    {
      buffer << mean(meanNeutSFS[i]) << '\t' << mean(meanRiskSFS[i]) << '\n';
    }
  gzFile gzout = gzopen(pars.outfile.c_str(),"w");
  gzwrite(gzout,buffer.str().c_str(),buffer.str().size());
  gzclose(gzout);
}

params parse_argv(int argc, char ** argv)
{
  using namespace boost::program_options;
  params rv;
  options_description desc("Calculate mean site frequency spectrum for non-risk and risk mutations separately");
  desc.add_options()
    ("help,h", "Produce help message")
    ("infile,i",value<string>(&rv.infile)->default_value(string()),"Input file name (binary, gzipped archive of simulated populations)")
    ("outfile,o",value<string>(&rv.outfile)->default_value(string()),"Output file name (will be gzipped)")
    ("nsam,n",value<unsigned>(&rv.nsam)->default_value(100u),"sample size (2x the number of diploids to sample")
    ("seed,S",value<unsigned>(&rv.seed)->default_value(0u),"Random number seed")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if(argc == 1 || vm.count("help"))
    {
      cerr << desc << '\n';
      exit(EXIT_SUCCESS);
    }

  if(! vm.count("infile") )
    {
      cerr << "Erorr: infile (--infile-i) not specified\n";
      exit(EXIT_FAILURE);
    }

  if(! vm.count("outfile") )
    {
      cerr << "Erorr: outfile (--outfile-i) not specified\n";
      exit(EXIT_FAILURE);
    }

  return rv;
}
