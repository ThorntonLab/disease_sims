/*
  Read phenotypes file and pop file.

  Outputs:
  1. Genetic component of pheno
  2. Random component of pheno
  3. # deleterious mutations on each haplotype
  4. The effect sizes of those mutations
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>

#include <fwdpp/diploid.hh>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>

#include <mutation_with_age.hpp>
#include <locking_routines.hpp>

using namespace std;
using namespace KTfwd;
using namespace boost::program_options;

struct params
{
  string ofile;
  vector<string> popfiles;
  unsigned nsam,seed;
};

params process_command_line(int argc, char ** argv);

int main( int argc, char ** argv )
{
  params options = process_command_line(argc,argv);  

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,options.seed);

  vector< unsigned > sfs_n(options.nsam-1,0u),sfs_c(options.nsam-1,0u);
  unsigned rep = 0;
  for( unsigned p = 0 ;p < options.popfiles.size() ; ++p )
    {
      ifstream popstream( options.popfiles[p].c_str() );
      unsigned nreps = 0;
      popstream.seekg(0,ios::end);
      long eostream = popstream.tellg();
      popstream.seekg(0,ios::beg);
      do
	{
	  glist gametes;
	  mlist mutations;
	  vector< pair< glist::iterator,glist::iterator > > diploids;
	  read_binary_pop( &gametes, &mutations, &diploids, boost::bind(mreader(),_1),popstream );
	  ++rep;

	  //Take a sample
	  pair<vector< pair<double,string> >,
	    vector< pair<double,string> > > sample = ms_sample_separate(r,&diploids,options.nsam);
	  
	  for( unsigned i = 0 ; i < sample.first.size() ; ++i )
	    {
	      sfs_n[ count( sample.first[i].second.begin(),
			    sample.first[i].second.end(),'1' ) - 1 ]++;
	    }
	  for( unsigned i = 0 ; i < sample.second.size() ; ++i )
	    {
	      sfs_c[ count( sample.second[i].second.begin(),
			    sample.second[i].second.end(),'1' ) - 1 ]++;
	    }
	}
      while( popstream.tellg() < eostream );
      popstream.close();
    }

  ostringstream buffer;
  for( unsigned i = 0 ; i < options.nsam-1 ; ++i )
    {
      buffer << (i+1) << '\t'
	     << double(sfs_n[i])/double(rep) << '\t'
	     << double(sfs_c[i])/double(rep) << '\n';
    }

  if ( options.ofile.empty() )
    {
      cout << buffer.str();
    }
  else
    {
      ofstream out(options.ofile.c_str());
      if(!out)
	{
	  cerr << "Error, " << options.ofile << " could not be opened for writing\n";
	  exit(10);
	}
      out << buffer.str();
      out.close();
    }
  exit(0);
}

params process_command_line(int argc, char ** argv)
{
  params rv;

  options_description desc("Process output from TFL2013_ind and generate case/control panel");

  desc.add_options()
    ("help,h", "Produce help message")
    ("popfile,p",value<vector<string> >(&rv.popfiles),"Population file output by TFL2013_ind.  May be used multiple times")
    ("nsam,n",value<unsigned>(&rv.nsam)->default_value(50),"Sample size")
    ("ofile,o",value<string>(&rv.ofile)->default_value(string()),"Name of output file")
    ("seed,S",value<unsigned>(&rv.seed)->default_value(0u),"random number seed")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if(argc == 1 || vm.count("help"))
    {
      cerr << desc << '\n';
      exit(0);
    }

  return rv;
}
