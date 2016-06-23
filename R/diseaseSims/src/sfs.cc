#include <Rcpp.h>
#include <zlib.h>
#include <functional>
#include <random>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <map>
#include <diseaseSims/mutation_with_age.hpp>
#include <diseaseSims/util.hpp>
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

  unsigned long NROW = 0;

  vector<unsigned> allID,allC,allnSFS,allrSFS;
  do 
    {
      //mlist mutations;
      //glist gametes;
      //dipvector diploids;
      //KTfwd::read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
      auto pop = readPop(gzin);
      unsigned N = (n < 2*pop.diploids.size()) ? n : 2*pop.diploids.size();
      NROW += N;
      vector<unsigned> rd(pop.diploids.size());
      unsigned dummy = 0;
      for_each(rd.begin(),rd.end(),[&dummy](unsigned & __u){__u=dummy++;});
      shuffle(rd.begin(),rd.end(),generator);
      //vector< pair<const mlist::iterator,unsigned> >  sampleRiskCounts,sampleNeutralCounts;
      map<unsigned,unsigned> sampleRiskCounts,sampleNeutralCounts;
      for( unsigned i = 0 ; i < N/2 ; ++i )
	{
	  assert( rd[i] < diploids.size() );
	  for(auto && m : pop.gametes[pop.diploids[rd[i]].first].mutations) sampleNeutralCounts[m]++;
	  for(auto && m : pop.gametes[pop.diploids[rd[i]].second].mutations) sampleNeutralCounts[m]++;
	  for(auto && m : pop.gametes[pop.diploids[rd[i]].first].smutations) sampleRiskCounts[m]++;
	  for(auto && m : pop.gametes[pop.diploids[rd[i]].second].smutations) sampleRiskCounts[m]++;
	}

      std::vector<unsigned> ID(N,REP),C(N),nSFS(N,0.),rSFS(N,0.);
      dummy=1;
      for_each(C.begin(),C.end(),[&dummy](unsigned & __u){__u=dummy++;});
      //for( unsigned i = 0 ; i < sampleNeutralCounts.size() ; ++i )
      for( auto && i : sampleNeutralCounts )
	{
	  nSFS[ i.second -1 ] += i.first;
	}
      //for( unsigned i = 0 ; i < sampleRiskCounts.size() ; ++i )
      for(auto && i : sampleRiskCounts )
	{
	  rSFS[ i.second - 1 ] += i.first;
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
