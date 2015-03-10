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
#include <gsl/gsl_randist.h>

using namespace std;
using namespace boost::accumulators;

using mean_acc = accumulator_set<double, stats<tag::mean> >;

double getVG( gzFile gzpin ); //not used...

int main( int argc, char ** argv)
{
  int argn=1;
  const char * popfilename = argv[argn++];
  const unsigned nreps = atoi(argv[argn++]);

  gzFile gzin = gzopen(popfilename,"rb");

  map<unsigned, pair<unsigned,mean_acc> > data;
  unsigned nrisk = 0;
  unsigned popsize= 0;
  for( unsigned i = 0 ; i < nreps ; ++i )
    {
      mlist mutations;
      glist gametes;
      dipvector diploids;
      KTfwd::read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
      popsize = diploids.size();
      for( auto mitr = mutations.begin() ; mitr != mutations.end() ; ++mitr )
	{
	  if( ! mitr->neutral )
	    {
	      ++nrisk;
	      auto __i = data.find( mitr->n );
	      if(__i == data.end())
		{
		  data[mitr->n] = make_pair(1u,mean_acc());
		  data[mitr->n].second(pow(mitr->s,2.));
		}
	      else
		{
		  __i->second.second( pow( mitr->s,2. ) );
		  __i->second.first++;
		}
	    }
	}
    }
  double SUM = 0.;
  for( const auto & d : data )
    {
      double p = double(d.first)/(2.*double(popsize));
      double mv = mean(d.second.second);
      double xx =  mv*(double(d.second.first)/double(nrisk))*p*(1.-p);
      SUM += 0.5*xx;
      cout << p << ' ' << SUM << '\n';
    }
}

double getVG( gzFile gzpin )
{
  unsigned nvals;
  gzread(gzpin,&nvals,sizeof(unsigned));
  double val;
  accumulator_set<double, stats<tag::variance> > av;
  for( unsigned i = 0 ; i < 2*nvals ; ++i )
    {
      gzread(gzpin,&val,sizeof(double));
      if( i%2==0. ) av(val);
    }
  return variance(av);
}
