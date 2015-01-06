//Site freq spectrum from effects file

#include <string>
#include <fstream>
#include <sstream>
#include <map>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>

#include <readSimOutput.hpp>

using namespace std;
using namespace boost::program_options;

typedef map<double,double> sfs;


struct params
{
  string indexfile,effectfile,outfile;
  bool normalize;
};

params process_argv( int argc, char ** argv );

void update_detail( const effectFileData & e,
		    sfs & meansfs,
		    const bool & normalize, 
		    const size_t & nsites)
{
  sfs::iterator itr = meansfs.find( e.count );
  if (itr == meansfs.end())
    {
      meansfs[e.count]= (normalize) ? 1./double(nsites) : 1.;
    }
  else
    {
      itr->second += (normalize) ? 1./double(nsites) : 1.;
    }
}

void update( const effectFileData & e,
	     sfs & mean_neut, 
	     sfs & mean_caus,
	     const bool & normalize,
	     const size_t & nsites)
{
  if ( e.esize != 0. )
    {
      update_detail( e, mean_caus , normalize, nsites );
    }
  else
    {
      update_detail( e, mean_neut , normalize, nsites );
    }
}

int main( int argc, char ** argv )
{
  params options = process_argv(argc,argv);

  cerr << options.normalize << '\n';
  gzFile effectsin = gzopen(options.effectfile.c_str(),"r");
  // ifstream estream( options.effectfile.c_str() );
  // if( ! estream )
  //   {
  //     cerr << "Error, " << options.effectfile << " could not be opened for reading\n";
  //     exit(10);
  //   }

  //data structures
  sfs meansfs_neut,meansfs_caus;

  unsigned nreps = 0;
  //estream.seekg(0,ios::end);
  //long eostream = estream.tellg();
  //estream.seekg(0,ios::beg);
  do
    {
      vector<effectFileData> effects = read_effect_file(effectsin);
      ++nreps;

      //sfs sfs_i;
      for( unsigned i = 0 ; i < effects.size() ; ++i )
	{
	  update( effects[i],meansfs_neut,meansfs_caus,options.normalize,effects.size() );
	}
    }
  while(!gzeof(effectsin));
  gzclose(effectsin);
  //while( estream.tellg() < eostream );

  ostringstream obuffer;
  
  //A little ugly.
  set<unsigned> ucounts;
  for( sfs::const_iterator i = meansfs_neut.begin() ; i != meansfs_neut.end() ;++i )
    {
      ucounts.insert( i->first );
    }
  for( sfs::const_iterator i = meansfs_caus.begin() ; i != meansfs_caus.end() ;++i )
    {
      ucounts.insert( i->first );
    }
  
  for( set<unsigned>::const_iterator i = ucounts.begin() ; i != ucounts.end() ; ++i )
    {
      obuffer << *i << '\t';
      sfs::const_iterator itr = meansfs_neut.find(*i);
      if( itr == meansfs_neut.end() )
	{
	  obuffer << "0\t";
	}
      else
	{
	  obuffer << itr->second/double(nreps) << '\t';
	}
      itr = meansfs_caus.find(*i);
      if( itr == meansfs_caus.end() )
	{
	  obuffer << "0\n";
	}
      else
	{
	  obuffer << itr->second/double(nreps) << '\n';
	}
    }

  if( options.outfile.empty() )
    {
      cout << obuffer.str();
    }

  exit(0);
}

params process_argv( int argc, char ** argv )
{
  params rv;

  options_description desc("Output mean sfs from effect size file.");

  desc.add_options()
    ("index,i",value<string>(&rv.indexfile)->default_value(string()),"Index file for simulation records")
    ("esizes,e",value<string>(&rv.effectfile)->default_value(string()),"Effects file from simulation")
    ("ofile,o",value<string>(&rv.outfile)->default_value(string()),"Output file name")
    ("normalize,n","Normalizes sfs by total no. of segregating mutations")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if(argc == 1 || vm.count("help"))
    {
      cerr << desc << '\n';
      exit(0);
    }

  if( vm.count("normalize") )
    {
      rv.normalize = true;
    }
  else
    {
      rv.normalize = false;
    }
       
  return rv;
}

