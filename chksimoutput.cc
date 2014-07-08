//Validates simulation output using the index file
//So far only checks number of mutations

#include <mutation_with_age.hpp>
//#include <TFL_fitness_models.hpp>
#include <simindex.hpp>
#include <fwdpp/IO.hpp>
#include <zlib.h>
#include <boost/bind.hpp>
#include <iostream>
using namespace std;

int main( int argc, char ** argv )
{
  int argn=1;
  const char * indexfile = argv[argn++];
  const char * popfile = argv[argn++];
  const char * effectfile = argv[argn++];

  simindex idx(indexfile);

  gzFile gzin=gzopen(popfile,"r"),gzin2=gzopen(effectfile,"r");

  for( unsigned i = 1 ; i <= idx.size() ; ++i )
    {
      mlist mutations;
      glist gametes;
      boost::container::vector< pair<glist::iterator,glist::iterator> > diploids;
      cerr << "checking record "
	  << i 
	  << " at positions "
	  << idx.hoffset(i) 
	  << " and "
	  << idx.eoffset(i) 
	  << '\n';
      gzseek( gzin, idx.hoffset(i), 0 );
      KTfwd::read_binary_pop( &gametes, &mutations, &diploids, boost::bind(gzmreader(),_1),gzin );

      gzseek( gzin2, idx.eoffset(i), 0 );
      unsigned nmuts;
      gzread( gzin2, &nmuts, sizeof(unsigned) );
      double evals[4];
      std::vector<double> pos;
      for( unsigned e = 0 ; e < nmuts ; ++e )
	{
	  gzread(gzin2,&evals[0],4*sizeof(double));
	  pos.push_back( evals[0] );
	}

      bool FAIL = false;
      if( pos.size() != mutations.size() )
	{
	  cerr << "Record "
	       << i 
	       << " mutation numbers don't match up: " << pos.size() << ' ' << mutations.size() << '\n';
	  FAIL = true;
	}
      else
	{
	  for( mlist::const_iterator itr = mutations.begin() ; 
	       itr != mutations.end() ; ++itr )
	    {
	      std::vector<double>::const_iterator pitr = find(pos.begin(),
							      pos.end(),
							      itr->pos);
	      if( pitr == pos.end() )
		{
		  cerr << "Record "
		       << i
		       << " mutation at position " << itr->pos 
		       << " in population is not found in effects file!\n";
		  FAIL = true;
		}
	    }
	}
      if(! FAIL )
	{
	  cerr << "Record "
	       << i 
	       << " checks out ok\n";
	}
    }
}
