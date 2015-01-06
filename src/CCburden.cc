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
#include <fwdpp/fwd_functional.hpp>

//#include <boost/bind.hpp>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/program_options.hpp>

#include <readCC.hpp>
//#include <locking_routines.hpp>

using namespace std;
using namespace KTfwd;
using namespace boost::interprocess;
using namespace boost::program_options;

struct params
{
  string ccfile,indexfile,ofile,esizefile,esizefileOUT,simindex;
  unsigned record_no,cases,controls;
  bool header;
};

params process_command_line(int argc, char ** argv);

struct minfo
{
  double s;
  double count,age;
};


int main( int argc, char ** argv )
{
  params options = process_command_line(argc,argv);  

  //Read effects
  ifstream index(options.simindex.c_str());
  if(!index)
    {
      cerr << "Error, cannot open " << options.simindex << " for reading\n";
      exit(10);
    }
  long pop_offset,pheno_offset,effect_offset;
  unsigned ith_rep;
  bool found = false;
  //get offset of the population and the phenotypes from the indexfile
  while( !found && !index.eof() )
    {
      index >> ith_rep >> effect_offset >> pheno_offset >> pop_offset >> ws;
      if( ith_rep == options.record_no )
	{
	  found = true;
	  break;
	}
    }
  ifstream ein(options.esizefile.c_str());
  if (! ein )
    {
      cerr << "Error, cannot open " << options.esizefile << " for reading\n";
      exit(10);
    }
  ein.seekg(effect_offset);
  unsigned nmuts;
  ein.read( reinterpret_cast<char *>(&nmuts),sizeof(unsigned) );
  map<double,minfo> mutinfo;
  for( unsigned i = 0 ; i < nmuts ; ++i )
    {
      double pos;
      minfo mi;
      ein.read( reinterpret_cast<char *>(&pos),sizeof(double) );
      ein.read( reinterpret_cast<char *>(&mi.s),sizeof(double) );
      ein.read( reinterpret_cast<char *>(&mi.count),sizeof(double) );
      ein.read( reinterpret_cast<char *>(&mi.age),sizeof(double) );
      mutinfo[pos]=mi;
    }
  ein.close();
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
  ostringstream effectbuffer;
  for( unsigned i = 0 ; i < options.controls ; ++i )
    {
      burdenbuffer << options.record_no << '\t'
		   << 0 << '\t'
		   << data.phenos[i].first << '\t'
		   << data.phenos[i].second << '\t'
		   << data.burden[i].first << '\t'
		   << data.burden[i].second << '\n';
      //mutation effect info
      for(unsigned mut = data.SN ; mut < data.SN + data.SC ; ++mut )
	{
	  if( data.geno_matrix[i][mut]==1 ) //hap 1
	    {
	      double pos = data.caus_pos[mut - data.SN];
	      map<double,minfo>::const_iterator mitr = mutinfo.find(pos);
	      assert( mitr != mutinfo.end() );
	      effectbuffer << options.record_no << '\t'
			   << i << '\t'
			   << "co" << '\t'
			   << mitr->second.s << '\t'
			   << mitr->second.count << '\t'
			   << mitr->second.age << '\n';
	    }
	  /*
	  if( data.geno_matrix[2*i + 1][mut]==1 ) //hap 2
	    {
	      double pos = data.caus_pos[mut - data.SN];
	      vector<TFLmtype>::const_iterator mitr = find_if(mvector.begin(),
							      mvector.end(),
							      std::bind(KTfwd::mutation_at_pos(),std::placeholders::_1,pos));
	      effectbuffer << i << '\t'
			   << "co" << '\t'
			   << 1 << 't'
			   << mitr->s << '\t'
			   << mitr->n << '\n';
	    }
	  */
	}
    }

  for( unsigned i = options.controls ; i < options.cases+options.controls ; ++i )
    {
      burdenbuffer << options.record_no << '\t'
		   << 1 << '\t'
		   << data.phenos[i].first << '\t'
		   << data.phenos[i].second << '\t'
		   << data.burden[i].first << '\t'
		   << data.burden[i].second << '\n';
      for(unsigned mut = data.SN ; mut < data.SN + data.SC ; ++mut )
	{
	  if( data.geno_matrix[i][mut]==1 ) 
	    {
	      double pos = data.caus_pos[mut - data.SN];
	      map<double,minfo>::const_iterator mitr = mutinfo.find(pos);
	      assert( mitr != mutinfo.end() );
	      effectbuffer << options.record_no << '\t'
			   << i << '\t'
			   << "ca" << '\t'
			   << mitr->second.s << '\t'
			   << mitr->second.count << '\t'
			   << mitr->second.age << '\n';
	    }
	  /*
	  if( data.geno_matrix[2*i + 1][mut]==1 ) //hap 2
	    {
	      double pos = data.caus_pos[mut - data.SN];
	      vector<TFLmtype>::const_iterator mitr = find_if(mvector.begin(),
							      mvector.end(),
							      std::bind(KTfwd::mutation_at_pos(),std::placeholders::_1,pos));
	      effectbuffer << i << '\t'
			   << "co" << '\t'
			   << 1 << 't'
			   << mitr->s << '\t'
			   << mitr->n << '\n';
	    }
	  */
	}
    }
  
  //Write output 
  if ( ! options.ofile.empty())
    {
      FILE * o_fh = fopen(options.ofile.c_str(),"a");
      int o_fd = fileno(o_fh);
      
      //flock o_lock = get_whole_flock();
      file_lock o_lock(options.ofile.c_str());
      //make sure our locking functions work..
      /*.
      assert( o_lock.l_type == F_WRLCK );
      assert( o_lock.l_whence == SEEK_SET );
      assert( o_lock.l_start == 0 );
      assert( o_lock.l_len == 0 );
      if (fcntl(o_fd,F_SETLKW,&o_lock) == -1)
	{
	  cerr << "ERROR: could not obtain lock on " << options.ofile << '\n';
	  exit(10);
	  }
      */
      if( ::write(o_fd,burdenbuffer.str().c_str(),burdenbuffer.str().size()) == -1 )
	{
	  cerr << "Error writing buffer to " << options.ofile << '\n';
	  exit(errno);
	}

      ofstream eout( options.esizefileOUT.c_str(),ios_base::app );
      eout << effectbuffer.str();
      eout.close();

      /*
	o_lock.l_type = F_UNLCK;
	if (fcntl(o_fd,F_UNLCK,&o_lock) == -1)
	{
	cerr << "ERROR: could not relesae lock on " << options.ofile << '\n';
	exit(10);
	}
      */
      o_lock.unlock();
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
    ("simindx,I",value<string>(&rv.simindex)->default_value(string()),"Index file for simulation data")
    ("ccfile,p",value<string>(&rv.ccfile)->default_value(string()),"Case control data files")
    ("efile,e",value<string>(&rv.esizefile)->default_value(string()),"Effect size file for input")
    ("efileout,E",value<string>(&rv.esizefileOUT)->default_value(string()),"Effect size file for output")
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
