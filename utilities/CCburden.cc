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

#include <fwdpp/IO.hpp>
#include <fwdpp/sampling_functions.hpp>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>

#include <readCC.hpp>
#include <locking_routines.hpp>

using namespace std;
using namespace KTfwd;
using namespace boost::program_options;

struct params
{
  string ccfile,indexfile,ofile;
  unsigned record_no,cases,controls;
  bool header;
};

params process_command_line(int argc, char ** argv);

int main( int argc, char ** argv )
{
  params options = process_command_line(argc,argv);  

  //Read in the case/control data
  bool fail;
  CCblock data = read_CC_record(options.indexfile.c_str(),
				options.ccfile.c_str(),
				options.record_no,
				&fail);

  
  if ( fail )
    {
      cerr << "Error: could not read record " 
	   << options.record_no
	   << " from " << options.ccfile 
	   << " using index file " << options.indexfile 
	   << '\n';
    }
  //output some stuff that may be useful
  ostringstream burdenbuffer;
  if ( options.header )
    {
      burdenbuffer << "record\tcase\tG\tE\th1\th2\n";
    }
  for( unsigned i = 0 ; i < options.controls ; ++i )
    {
      burdenbuffer << options.record_no << '\t'
		   << 0 << '\t'
		   << data.phenos[i].first << '\t'
		   << data.phenos[i].second << '\t'
		   << data.burden[i].first << '\t'
		   << data.burden[i].second << '\n';
    }
  for( unsigned i = options.controls ; i < options.cases+options.controls ; ++i )
    {
      burdenbuffer << options.record_no << '\t'
		   << 1 << '\t'
		   << data.phenos[i].first << '\t'
		   << data.phenos[i].second << '\t'
		   << data.burden[i].first << '\t'
		   << data.burden[i].second << '\n';
    }

  //Write output 
  if ( ! options.ofile.empty())
    {
      FILE * o_fh = fopen(options.ofile.c_str(),"a");
      int o_fd = fileno(o_fh);
      
      flock o_lock = get_whole_flock();
      
      //make sure our locking functions work...
      assert( o_lock.l_type == F_WRLCK );
      assert( o_lock.l_whence == SEEK_SET );
      assert( o_lock.l_start == 0 );
      assert( o_lock.l_len == 0 );
      if (fcntl(o_fd,F_SETLKW,&o_lock) == -1)
	{
	  cerr << "ERROR: could not obtain lock on " << options.ofile << '\n';
	  exit(10);
	}
      if( ::write(o_fd,burdenbuffer.str().c_str(),burdenbuffer.str().size()) == -1 )
	{
	  cerr << "Error writing buffer to " << options.ofile << '\n';
	  exit(errno);
	}

      o_lock.l_type = F_UNLCK;
      if (fcntl(o_fd,F_UNLCK,&o_lock) == -1)
	{
	  cerr << "ERROR: could not relesae lock on " << options.ofile << '\n';
	  exit(10);
	}
      fflush(o_fh);
      fclose(o_fh);
    }
  else
    {
      cout << burdenbuffer.str();
    }
  exit(0);
}

params process_command_line(int argc, char ** argv)
{
  params rv;

  options_description desc("Process output from TFL2013_ind and generate case/control panel");

  desc.add_options()
    ("help,h", "Produce help message")
    ("indexfile,i",value<string>(&rv.indexfile)->default_value(string()),"Index file for case/control data")
    ("ccfile,p",value<string>(&rv.ccfile)->default_value(string()),"Case control data files")
    ("recordno,r",value<unsigned>(&rv.record_no)->default_value(numeric_limits<unsigned>::max()),"Record value to look up from input files")
    ("ofile,o",value<string>(&rv.ofile)->default_value(string()),"Name of output file")
    ("noheader,N","If used, no header information will be output")
    ("cases,C",value<unsigned>(&rv.cases)->default_value(0u),"Number of cases.")
    ("controls,c",value<unsigned>(&rv.controls)->default_value(0u),"Number of controls.")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if (vm.count("noheader"))
    {
      rv.header=false;
    }
  else
    {
      rv.header=true;
    }
  if(argc == 1 || vm.count("help"))
    {
      cerr << desc << '\n';
      exit(0);
    }

  return rv;
}
