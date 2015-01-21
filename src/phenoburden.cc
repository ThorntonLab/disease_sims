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

//#include <boost/bind.hpp>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/program_options.hpp>

#include <diseaseSims/mutation_with_age.hpp>
//#include <locking_routines.hpp>

using namespace std;
using namespace KTfwd;
using namespace boost::interprocess;
using namespace boost::program_options;

struct params
{
  string popfile,phenofile,indexfile,ofile;
  unsigned record_no;
  bool header;
};

params process_command_line(int argc, char ** argv);

int main( int argc, char ** argv )
{
  params options = process_command_line(argc,argv);  

  //lookup the record
  long pop_offset,pheno_offset,effect_offset;
  unsigned ith_rep;
  bool found = false;
  //get offset of the population and the phenotypes from the indexfile
  ifstream index( options.indexfile.c_str() );
  while( !found && !index.eof() )
    {
      index >> ith_rep >> effect_offset >> pheno_offset >> pop_offset >> ws;
      if( ith_rep == options.record_no )
	{
	  found = true;
	  break;
	}
    }

  if ( ! found )
    {
      cerr << "Error: record number " << options.record_no << " not found in " << options.indexfile << '\n';
      exit(10);
    }
  index.close();

  //Read in the population
  glist gametes;
  mlist mutations;
  vector< pair< glist::iterator,glist::iterator > > diploids;

  ifstream popstream( options.popfile.c_str() );
  popstream.seekg( pop_offset );
  read_binary_pop( &gametes, &mutations, &diploids, std::bind(mreader(),std::placeholders::_1),popstream );
  popstream.close();

  //read in the phenotypes
  ifstream phenostream(options.phenofile.c_str());
  phenostream.seekg(pheno_offset);
  vector< std::pair<double,double> > phenotypes;
  unsigned n;
  phenostream.read( reinterpret_cast<char *>(&n),sizeof(unsigned) );

  assert(n==diploids.size());
  double p[2];
  for(unsigned i=0;i<n;++i)
    {
      phenostream.read( reinterpret_cast<char *>(&p[0]),2*sizeof(double) );
      phenotypes.push_back( std::make_pair(p[0],p[1]) );
    }
  phenostream.close();

  //output some stuff that may be useful
  ostringstream burdenbuffer;
  if ( options.header )
    {
      burdenbuffer << "record diploid G E haplo esize count\n";
    }
  for( unsigned i = 0 ; i < diploids.size() ; ++i )
    {
      if( diploids[i].first->smutations.empty() )
	{
	  burdenbuffer << options.record_no << '\t'
		       << i << '\t' 
		       << phenotypes[i].first 
		       << '\t' << phenotypes[i].second << '\t'
		       << 0 << '\t' 
		       << "NA\tNA\n";
	}
      else
	{
	  for( gtype::mcont_const_iterator mitr = diploids[i].first->smutations.begin() ; 
	       mitr != diploids[i].first->smutations.end() ; ++mitr )
	    {
	      burdenbuffer << options.record_no << '\t'
			   << i << '\t' 
			   << phenotypes[i].first << '\t' 
			   << phenotypes[i].second << '\t'
			   << 0 << '\t' 
			   << (*mitr)->s << '\t'
			   << (*mitr)->n << '\n';
	    }
	}

      if( diploids[i].second->smutations.empty() )
	{
	  burdenbuffer << options.record_no << '\t'
		       << i << '\t'
		       << phenotypes[i].first << '\t' 
		       << phenotypes[i].second << '\t' 
		       << 1 << '\t' 
		       << "NA\tNA\n";
	}
      else
	{
	  for( gtype::mcont_const_iterator mitr = diploids[i].second->smutations.begin() ; 
	       mitr != diploids[i].second->smutations.end() ; ++mitr )
	    {
	      burdenbuffer << options.record_no << '\t'
			   << i << '\t' 
			   << phenotypes[i].first << '\t' 
			   << phenotypes[i].second << '\t'
			   << 1 << '\t' 
			   << (*mitr)->s << '\t'
			   << (*mitr)->n << '\n';
	    }
	}
    }

  //Write output 
  if ( ! options.ofile.empty())
    {
      FILE * o_fh = fopen(options.ofile.c_str(),"a");
      int o_fd = fileno(o_fh);
      
      file_lock o_lock(options.ofile.c_str());
      scoped_lock<file_lock> s_lock(o_lock);
      //flock o_lock = get_whole_flock();
      
      //make sure our locking functions work...
      /*
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

      /*
      o_lock.l_type = F_UNLCK;
      if (fcntl(o_fd,F_UNLCK,&o_lock) == -1)
	{
	  cerr << "ERROR: could not relesae lock on " << options.ofile << '\n';
	  exit(10);
	}
      */
      fflush(o_fh);
      o_lock.unlock();
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
    ("indexfile,i",value<string>(&rv.indexfile)->default_value(string()),"Index file output by TFL2013_ind")
    ("popfile,p",value<string>(&rv.popfile)->default_value(string()),"Population file output by TFL2013_ind")
    ("phenofile,P",value<string>(&rv.phenofile)->default_value(string()),"Phenotypes file output by TFL2013_ind")
    ("recordno,r",value<unsigned>(&rv.record_no)->default_value(numeric_limits<unsigned>::max()),"Record value to look up from input files")
    ("ofile,o",value<string>(&rv.ofile)->default_value(string()),"Name of output file")
    ("noheader,N","If used, no header information will be output")
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
