//Validates simulation output using the index file
//So far only checks number of mutations

#include <mutation_with_age.hpp>
#include <simindex.hpp>
#include <fwdpp/IO.hpp>
#include <zlib.h>
#include <boost/bind.hpp>
#include <iostream>
#include <string>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;

struct options
{
  string ifile,pfile,efile;
};

options pargs(int argc, char ** argv);

int main( int argc, char ** argv )
{
  int argn=1;
  options o = pargs(argc,argv);

  simindex idx(o.ifile.c_str());

  gzFile gzin=gzopen(o.pfile.c_str(),"r"),gzin2=gzopen(o.efile.c_str(),"r");

  for( unsigned i = 1 ; i <= idx.size() ; ++i )
    {
      mlist mutations;
      glist gametes;
      std::vector< pair<glist::iterator,glist::iterator> > diploids;
      cerr << "checking record "
	  << i 
	  << " at positions "
	  << idx.hoffset(i) 
	  << " and "
	  << idx.eoffset(i) 
	  << '\n';
      gzseek( gzin, idx.hoffset(i), 0 );
      KTfwd::read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );

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

options pargs(int argc, char ** argv)
{
  options rv;
  options_description desc("Reads through the simulation output and index file to validate the output.  Currently, the program reads the index file, the population file, and the effects file, and makes sure that all the expected mutations in the population file and effects file match up.");
  desc.add_options()
    ("help,h","Show help")
    ("index,i",value<string>(&rv.ifile),"Index file for simulation output")
    ("pop,p",value<string>(&rv.pfile),"Binary file containing population")
    ("effects,e",value<string>(&rv.efile),"Binary file containing mutation effects")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if(argc == 1 || vm.count("help") ||
     (!vm.count("index")||!vm.count("pop")||!vm.count("effects")))
    {
      cerr << desc << '\n';
      exit(0);
    }

  return rv;
}
