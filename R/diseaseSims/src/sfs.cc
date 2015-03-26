#include <Rcpp.h>
#include <zlib.h>
#include <functional>
#include <random>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>

#include <fwdpp/IO.hpp>
#include <diseaseSims/mutation_with_age.hpp>

using namespace std;

//' Site frequency spectra
//' @param popfile The file containing the simulated populations
//' @param n The desired sample size.  If n == 2N, the entire population is sampled
//' @param seed A random number seed
//' @return A data frame with 4 columns:
//' @return replicate = an arbitrary sample ID.  This is NOT the same ID number as the user may have assigned when running the simulation!!
//' @return i = the number of occurences of a mutation in the population/sample
//' @return n = the SFS for neutral variants
//' @return r = the SFS for risk variants
//[[Rcpp::export]]
Rcpp::DataFrame sfs( const std::string & popfile,
			 const unsigned & n,
			 const unsigned & seed = 0 )
{
  Rcpp::NumericMatrix rv;
  gzFile gzin = gzopen(popfile.c_str(),"rb");
  mt19937 generator;
  generator.seed(seed);
  unsigned REP = 0;

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
  unsigned long NROW = 0;

  vector<unsigned> allID,allC,allnSFS,allrSFS;
  do 
    {
      mlist mutations;
      glist gametes;
      dipvector diploids;
      KTfwd::read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
      unsigned N = (n < 2*diploids.size()) ? n : 2*diploids.size();
      NROW += N;
      vector<unsigned> rd(diploids.size());
      unsigned dummy = 0;
      for_each(rd.begin(),rd.end(),[&dummy](unsigned & __u){__u=dummy++;});
      shuffle(rd.begin(),rd.end(),generator);
      vector< pair<const mlist::iterator,unsigned> >  sampleRiskCounts,sampleNeutralCounts;
      for( unsigned i = 0 ; i < N/2 ; ++i )
	{
	  assert( rd[i] < diploids.size() );
	  for_each( diploids[ rd[i] ].first->mutations.begin(),diploids[ rd[i] ].first->mutations.end(),std::bind(updater,std::placeholders::_1,ref(sampleNeutralCounts)) );
	  for_each( diploids[ rd[i] ].second->mutations.begin(),diploids[ rd[i] ].second->mutations.end(),std::bind(updater,std::placeholders::_1,ref(sampleNeutralCounts)) );
	  for_each( diploids[ rd[i] ].first->smutations.begin(),diploids[ rd[i] ].first->smutations.end(),std::bind(updater,std::placeholders::_1,ref(sampleRiskCounts)) );
	  for_each( diploids[ rd[i] ].second->smutations.begin(),diploids[ rd[i] ].second->smutations.end(),std::bind(updater,std::placeholders::_1,ref(sampleRiskCounts)) );
	}

      std::vector<unsigned> ID(N,REP),C(N),nSFS(N,0.),rSFS(N,0.);
      dummy=0;
      for_each(C.begin(),C.end(),[&dummy](unsigned & __u){__u=dummy++;});
      for( unsigned i = 0 ; i < sampleNeutralCounts.size() ; ++i )
	{
	  assert( sampleNeutralCounts[i].second > 0 );
	  assert( sampleNeutralCounts[i].second <= N );
	  nSFS[ sampleNeutralCounts[i].second - 1 ]++;
	}
      for( unsigned i = 0 ; i < sampleRiskCounts.size() ; ++i )
	{
	  assert( sampleRiskCounts[i].second > 0 );
	  assert( sampleRiskCounts[i].second <= N );
	  rSFS[ sampleRiskCounts[i].second - 1 ]++;
	}

      ++REP;

      //Update the return value
      copy(ID.begin(),ID.end(),back_inserter(allID));
      copy(C.begin(),C.end(),back_inserter(allC));
      copy(nSFS.begin(),nSFS.end(),back_inserter(allnSFS));
      copy(rSFS.begin(),rSFS.end(),back_inserter(allrSFS));

      //Check if at EOF
      int c = gzgetc(gzin);
      if(gzeof(gzin)) break;
      else gzungetc(c,gzin);
    }
  while(! gzeof(gzin) );

  return Rcpp::DataFrame::create( Rcpp::Named("replicate") = allID,
				  Rcpp::Named("i") = allC,
				  Rcpp::Named("n") = allnSFS,
				  Rcpp::Named("r") = allrSFS );
}
