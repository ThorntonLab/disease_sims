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

struct fdata 
{
  unsigned nm; //number of mutations at this frequency
  mean_acc ez,ezsq;
  fdata() : nm(0u),ez(mean_acc()),ezsq(mean_acc()) {}
};

int main( int argc, char ** argv)
{
  int argn=1;
  const char * popfilename = argv[argn++];
  const unsigned nreps = atoi(argv[argn++]);

  gzFile gzin = gzopen(popfilename,"rb");

  //map<unsigned, pair<unsigned,mean_acc> > data;
  map<unsigned,fdata> data;
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
		  data[mitr->n] = fdata();
		  data[mitr->n].nm++;
		  data[mitr->n].ez(mitr->s);
		  data[mitr->n].ezsq(pow(mitr->s,2.));
		}
	      else
		{
		  __i->second.nm++;
		  __i->second.ez(mitr->s);
		  __i->second.ezsq( pow(mitr->s,2.) );
		}
	    }
	}
    }
  double SUM = 0.,SUM2=0.;
  for( const auto & d : data )
    {
      double x = double(d.first)/(2.*double(popsize));
      double ez_x = mean(d.second.ezsq);
      double fx = double(d.second.nm)/double(nrisk);
      SUM += 0.5*ez_x*fx*x*(1.-x);

      //double x2 = 0.;
      //for( const auto & d2 : data ) {
	//double p2 = double(d2.first)/(2.*double(popsize));
	  //x2 += mean(d.second.second)*mean(d2.second.second)*double(d2.second.first)/double(nrisk)*((d.first == d2.first) ? p*p : 2.*p*p2);
      //}
      //SUM2 += 0.5*(double(d.second.first)/double(nrisk))*x2;
      cout << x << ' ' << SUM << '\n';
      //cout << p << ' ' << SUM << ' ' << SUM2 << '\n';
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
