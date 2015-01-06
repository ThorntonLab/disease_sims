/*
  Make case/control "panel" assuming a liability threshold on phenotype.
  In this case, individuals with phenotypes >= the x-th quantile
  of the population's phenotypic distribution are potential cases.

  There are two output files:

  1. anovafile_index:  this file records the record_no passed to the command line,
  and the offset (bytes) of where that record begins in anovafile

  2. anovafile: This fine is written in binary format.  It contains a 
  series of records of case/control panels.  Each record has the following format:
  Four integers: number of controls, number of cases, NN = number of neutral mutations, NC = number of causative mutations
  NN + NC doubles, which are the positions of the mutations (neutral then causative mutations, respectively).

  A genotype can take on 1 of 3 values, and is equal to the # of copies of the minor allele:
  0 = homozygote for major allele
  1 = heterozygote
  2 = homozygote for minor allele

  For each control and each case, the genotype info are recorded in the following format:
  1.  N1 = An unsigned integer representing the number of 1 (het) genotypes
  2.  N1 unsigned integers representing the indexes where the 1 values are located
  3.  N2 = An unsigned integer representing the number of 2 (MA hom) genotypes
  4.  N1 unsigned integers representing the indexes where the 2 values are located

  These indexes are stored in the same order as the positions.

  See single_marker_test.R that comes with this distribution for how to read these data files in 
  in R.

  After the genotypes, there are 2*(# controls + # cases) integers.  For each case, and then
  for each control, there is a pair of integers representing the # of causative mutations
  on the "paternal" and "maternal" haplotype, respectively.

  Finally, there are 2*(# controls + # cases) doubles. For each individual, each pair of 
  doubles represents the genetic and random component of phenotype, respectively.

  Comments: The index file seeks only to the very beginning of a record.  However,
  once you know the # of controls. cases, and mutations, you can seek around within
  a record easily using standard C or C++ functions as appropriate.
 */

#include <mutation_with_age.hpp>
#include <ccintermediate.hpp>
#include <simindex.hpp>
//#include <locking_routines.hpp>

#include <Sequence/SimData.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <vector>
#include <cassert>
#include <algorithm>
#include <cstdlib> 
#include <set>

#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/program_options.hpp>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <fwdpp/IO.hpp>
#include <fwdpp/sampling_functions.hpp>


using namespace std;
using namespace Sequence;
using namespace KTfwd;
using namespace boost::interprocess;
using namespace boost::program_options;

struct params
{
  string indexfile,popfile,phenofile,anovafile,anova_indexfile,idfile;
  unsigned twoN,record_no,ncases,ncontrols,seed;
  double case_proportion,crange;
  bool gzinput,gzoutput;
  params();
  bool files_undef( void ) const;
  void report_empty( std::ostream & out ) const;
  bool params_ok( void ) const;
};

params::params() : indexfile(string()),
		   popfile(string()),
		   phenofile(string()),
		   anovafile(string()),
		   anova_indexfile(string()),
		   idfile(string()),
		   twoN(0),
		   record_no(std::numeric_limits<unsigned>::max()),
		   ncases(0),
		   ncontrols(0),
		   seed(0),
		   case_proportion(1.),
		   crange(0.5),
		   gzinput(false),
		   gzoutput(false)
{
}

params process_command_line(int argc, char ** argv);

std::pair<double,double> phenosums(const vector<pair<double,double> > & phenos, const double & case_proportion, double * cutoff);

int main(int argc, char ** argv)
{
  params options = process_command_line(argc,argv);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,options.seed);
  
  //Read in the population
  glist gametes;
  mlist mutations;
  vector< pair< glist::iterator,glist::iterator > > diploids;

  bool found = false;
  //get offset of the population and the phenotypes from the indexfile
  simindex index(options.indexfile.c_str());

  if( index.file_problem() )
    {
      cerr << "Error reading index file\n";
      exit(10);
    }

  if( index.eexists( options.record_no ) )
    {
      found = true;
    }
 
  if ( ! found )
    {
      cerr << "Error: record number " << options.record_no << " not found in " << options.indexfile << '\n';
      exit(10);
    }

  if ( options.gzinput )
    {
      gzFile gzin = gzopen(options.popfile.c_str(),"rb");
      gzseek( gzin, index.hoffset(options.record_no), 0);
      read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
      gzclose(gzin);
    }
  else
    {
      ifstream popstream( options.popfile.c_str() );
      popstream.seekg( index.hoffset(options.record_no) );
      read_binary_pop( &gametes, &mutations, &diploids, std::bind(mreader(),std::placeholders::_1),popstream );
      popstream.close();
    }

  if ( gametes.size() > 2*diploids.size() )
    {
      cerr << "Error: there are " << gametes.size() 
	   << " gametes for only " << diploids.size() 
	   << " diploids.  Uncool!\n";
      exit(10);
    }
  //Make sure that ncases + ncontrols < N
  if ( (options.ncontrols + options.ncases) > diploids.size() )
    {
      std::cerr << "Error: population size is " << diploids.size() << " diploids. "
		<< "Sum of cases and controls is " << (options.ncontrols + options.ncases) << '\n';
      exit(10);
    }

  //Read in the phenotype data
  vector< pair<double,double> > phenotypes;

  if( options.gzinput )
    {
      gzFile gzin = gzopen( options.phenofile.c_str(),"rb" );
      gzseek( gzin, index.poffset(options.record_no), 0);
      unsigned nphenos;
      gzread(gzin,&nphenos,sizeof(unsigned));
      for(unsigned i = 0 ; i < nphenos ; ++i )
	{
	  //x is the genetic contribution to phenotype. y is the Gaussian noise from the simulation.
	  //Phenotype of the individual is x+y
	  double x,y;
	  gzread(gzin,&x,sizeof(double));
	  gzread(gzin,&y,sizeof(double)); 
	  phenotypes.push_back( make_pair(x,y) );
	}
      gzclose(gzin);
    }
  else
    {
      ifstream phenostream( options.phenofile.c_str() );
      phenostream.seekg( index.poffset(options.record_no) );
      
      unsigned nphenos;
      phenostream.read( reinterpret_cast<char *>(&nphenos), sizeof(unsigned) );
      for( unsigned i = 0 ; i < nphenos ; ++i )
	{
	  //x is the genetic contribution to phenotype. y is the Gaussian noise from the simulation.
	  //Phenotype of the individual is x+y
	  double x,y;
	  phenostream.read( reinterpret_cast<char *>(&x), sizeof(double) );
	  phenostream.read( reinterpret_cast<char *>(&y), sizeof(double) );
	  phenotypes.push_back( make_pair(x,y) );
	}
      phenostream.close();
    }

  assert( phenotypes.size() == diploids.size() );
  //The real work starts here

  //1. get mean, sd, and upper quantile of pheno distribution
  double cutoff;
  pair<double,double> mean_sd = phenosums(phenotypes,options.case_proportion,&cutoff);

  /*
    2. Assign an individual to be a putative case, control, nor not included
    
    putative case = phenotype >= cutoff
    putative control = phenotype within mean +/- sd of population distribution
  */
  vector< size_t > put_cases,put_controls;

  for( size_t i = 0 ; i < phenotypes.size() ; ++i )
    {
      const double P = phenotypes[i].first+phenotypes[i].second;
      if( P >= cutoff )
	{
	  put_cases.push_back(i);
	}
      /*
	Issue alert!!!
	TFL (2013) claim that "controls are within 1 standard
	deviation of the mean", by which they mean in terms of phenotype.
	
	Well, that was technically true, but a problem for reproducibility
	because our code to make case/control panels required that a putative
	control's phenotype be within 0.5*sd of population mean phenotype!!!
      */
      else if ( P >= mean_sd.first - options.crange*mean_sd.second &&
		P <= mean_sd.first + options.crange*mean_sd.second )
	{
	  put_controls.push_back(i);
	}
    }

  /*
    Check: do putative cases and controls overlap?

    If so, can we get rid of the overlap and still have sufficient #s in each category?

    If not, we warn that the CC data may be invalid, as cases and controls will share individuals.
  */
  vector< size_t > isect;
  sort( put_controls.begin(), put_controls.end() );
  sort( put_cases.begin(), put_cases.end() );
  set_intersection( put_controls.begin(),
		    put_controls.end(),
		    put_cases.begin(),
		    put_cases.end(),
		    std::back_inserter(isect) );
  if ( ! isect.empty() )
    {
      //Then case & control lists overlap and share individuals
      ostringstream wbuffer;
      wbuffer << "Warning: list of putative cases and controls have individuals in common.\n"; 
      //Get the set differences now
      vector< size_t > ucontrols, ucases;
      set_difference( put_controls.begin() , put_controls.end(),
		      put_cases.begin(), put_cases.end(),
		      std::back_inserter( ucontrols ) );
      set_difference( put_cases.begin() , put_cases.end(),
		      put_controls.begin(), put_controls.end(),
		      std::back_inserter( ucases ) );
      bool we_can_fix_this = true;
      if ( ucontrols.size() < options.ncontrols )
	{
	  we_can_fix_this = false;
	  wbuffer << "Warning: there are too few unique putative controls. "
		  << "There are " << options.ncontrols << " controls desired, but only " << ucontrols.size() << " "
		  << "diploids do not overlap with list of putative cases.\n";
	}
      if ( ucases.size() < options.ncases )
	{
	  we_can_fix_this = false;
	  wbuffer << "Warning: there are too few unique putative cases. "
		  << "There are " << options.ncases << " cases desired, but only " << ucases.size() << " "
		  << "diploids do not overlap with list of putative controls.\n";
	}
      if ( we_can_fix_this )
	{
	  put_cases = ucases;
	  put_controls = ucontrols;
	  ucases.clear();
	  ucontrols.clear();
	}
      else
	{
	  cerr << wbuffer.str()
	       << "Warning: your cases and controls may share individuals for record number " << options.record_no << "\n."
	       << "We are proceeding anyways...\n";
	}
    }

  //Randomize lists just for fun
  std::function< size_t (size_t) > rand = std::bind(&gsl_ran_flat, r, 0,double(put_cases.size()));
  random_shuffle(put_cases.begin(),put_cases.end(),rand);
  rand = std::bind(&gsl_ran_flat, r, 0,double(put_controls.size()));
  random_shuffle(put_controls.begin(),put_controls.end(),rand);

  cc_intermediate * ccblocks = new cc_intermediate(  process_population(diploids,phenotypes,
									put_controls,
									put_cases,
									options.ncontrols,
									options.ncases) );
  assert( ccblocks->min_n.size() == ccblocks->neutral.numsites() );
  assert( ccblocks->min_c.size() == ccblocks->causative.numsites() );
  
  //free up memory
  
  /*
  cerr << "deleting dips " << diploids.size() << ' ';
  //diploids.clear();
  cerr << "gams " << gametes.size() << ' ';
  delete gametes;
  //gametes.clear();
  cerr << "muts ";
  delete mutations;
  //mutations.clear();
  cerr << "done\n";
  */
  //3. Output to buffers

  //first, the ccdata to a buffer
  ostringstream ccbuffer;
  ccbuffer.write( reinterpret_cast<char *>(&options.ncontrols),sizeof(unsigned) );
  ccbuffer.write( reinterpret_cast<char *>(&options.ncases),sizeof(unsigned) );
  unsigned temp = ccblocks->neutral.numsites();
  ccbuffer.write( reinterpret_cast<char *>(&temp),sizeof(unsigned) );
  temp = ccblocks->causative.numsites();
  ccbuffer.write( reinterpret_cast<char *>(&temp),sizeof(unsigned) );

  for( SimData::const_pos_iterator p = ccblocks->neutral.pbegin() ; 
       p < ccblocks->neutral.pend() ; ++p )
    {
      double x = *p;
      ccbuffer.write( reinterpret_cast<char *>(&x),sizeof(double) );
    }

  for( SimData::const_pos_iterator p = ccblocks->causative.pbegin() ; 
       p < ccblocks->causative.pend() ; ++p )
    {
      double x = *p;
      ccbuffer.write( reinterpret_cast<char *>(&x),sizeof(double) );
    }

  //iterate over the diploids and write block for association tests
  for( unsigned ind = 0 ; ind < ccblocks->neutral.size() ; ind += 2 ) 
    {
      vector<unsigned> ones,twos;
      //neutral genotypes for this individual
      for( unsigned site = 0 ; site < ccblocks->neutral.numsites() ; ++site )
	{
	  //count # copies of minor allele at this site in this individual
	  unsigned cminor = ( (ccblocks->neutral[ind][site] == ccblocks->min_n[site]) ? 1 : 0 ) +
	    ( (ccblocks->neutral[ind+1][site] == ccblocks->min_n[site]) ? 1 : 0 );
	  if( cminor==1 )
	    {
	      ones.push_back(site);
	    }
	  else if ( cminor == 2 )
	    {
	      twos.push_back(site);
	    }
	  //ccbuffer.write( reinterpret_cast<char *>(&cminor), sizeof(unsigned) );
	}
      //causative genotypes for this individual
      for( unsigned site = 0 ; site < ccblocks->causative.numsites() ; ++site )
	{
	  //count # copies of minor allele at this site in this individual
	  unsigned cminor = ( (ccblocks->causative[ind][site] == ccblocks->min_c[site]) ? 1 : 0 ) +
	    ( (ccblocks->causative[ind+1][site] == ccblocks->min_c[site]) ? 1 : 0 );
	  //ccbuffer.write( reinterpret_cast<char *>(&cminor), sizeof(unsigned) );
	  if( cminor==1 )
	    {
	      ones.push_back(ccblocks->neutral.numsites()+site);
	    }
	  else if ( cminor == 2 )
	    {
	      twos.push_back(ccblocks->neutral.numsites()+site);
	    }
	}
      //update the buffer
      unsigned n = ones.size();
      ccbuffer.write( reinterpret_cast<char *>(&n),sizeof(unsigned) );
      ccbuffer.write( reinterpret_cast<char *>(&ones[0]),n*sizeof(unsigned) );
      n = twos.size();
      ccbuffer.write( reinterpret_cast<char *>(&n),sizeof(unsigned) );
      ccbuffer.write( reinterpret_cast<char *>(&twos[0]),n*sizeof(unsigned) );
    }

  //Now, output # of causative mutations on each haplotype carried by this diploid
  for( unsigned ind = 0 ; ind < ccblocks->causative.size() ; ind += 2 ) 
    {
      unsigned ncaus = count( ccblocks->causative[ind].begin(),
			      ccblocks->causative[ind].end(),'1' );
      ccbuffer.write( reinterpret_cast<char *>(&ncaus), sizeof(unsigned) );
      ncaus = count( ccblocks->causative[ind+1].begin(),
		     ccblocks->causative[ind+1].end(),'1' );
      ccbuffer.write( reinterpret_cast<char *>(&ncaus), sizeof(unsigned) );
    }
  //Give phenos of controls & cases
  for( vector< pair<double,double> >::const_iterator i = ccblocks->phenotypes.begin() ; 
       i != ccblocks->phenotypes.end() ; ++i )
    {
      double x = i->first;
      ccbuffer.write( reinterpret_cast<char *>(&x),sizeof(double) );
      x = i->second;
      ccbuffer.write( reinterpret_cast<char *>(&x),sizeof(double) );
    }

  //free up RAM
  delete ccblocks;

  //4. Output to files

  //obtain file lock on index ASAP
  FILE * ai_fh = fopen(options.anova_indexfile.c_str(),"a");
  int ai_fd = fileno(ai_fh);
  
  file_lock ai_lock(options.anova_indexfile.c_str());
  //struct flock ai_lock = get_whole_flock();
  
  //make sure our locking functions work...
  /*
  assert( ai_lock.l_type == F_WRLCK );
  assert( ai_lock.l_whence == SEEK_SET );
  assert( ai_lock.l_start == 0 );
  assert( ai_lock.l_len == 0 );
  if (fcntl(ai_fd,F_SETLKW,&ai_lock) == -1)
    {
      cerr << "ERROR: could not obtain lock on " << options.anova_indexfile << '\n';
      exit(10);
    }
  */
  ostringstream idbuffer;
  if( !options.idfile.empty() ) //then we want to write the individual id indexes
    {
      idbuffer.write( reinterpret_cast<char *>(&options.ncontrols), sizeof(unsigned) );  
      idbuffer.write( reinterpret_cast<char *>(&options.ncases), sizeof(unsigned) );  
      for(unsigned i = 0 ; i < options.ncontrols ; ++i )
	{
	  idbuffer.write( reinterpret_cast<char *>(&put_controls[i]), sizeof(unsigned) );  
	}
      for(unsigned i = 0 ; i < options.ncases ; ++i )
	{
	  idbuffer.write( reinterpret_cast<char *>(&put_cases[i]), sizeof(unsigned) );  
	}
    }

  if( options.gzoutput )
    {
      gzFile gzout = gzopen( options.anovafile.c_str(),"a" );
      int written = gzwrite( gzout, ccbuffer.str().c_str(), ccbuffer.str().size() );
      if(!written)
	{
	  cerr << "Error writing to " << options.anovafile << '\n';
	  exit(10);
	}
      //gzflush(gzout,Z_FINISH);
      gzclose(gzout);

      int idwritten = 0;
      if( !options.idfile.empty() ) //then we want to write the individual id indexes
	{
	  gzout = gzopen( options.idfile.c_str(), "a" );
	  idwritten = gzwrite( gzout, idbuffer.str().c_str(), ccbuffer.str().size() );
	  if(! idwritten)
	    {
	      cerr << "Error writing to " << options.idfile << '\n';
	      exit(10);
	    }
	  gzclose(gzout);
	}
      ostringstream buffer;
      buffer << options.record_no << ' ' << written;
      if( !options.idfile.empty() )
	{
	  buffer << ' ' << idwritten << '\n';
	}
      else 
	{
	  buffer << '\n';
	}
      if( ::write(ai_fd,buffer.str().c_str(),buffer.str().size()) == -1 )
	{
	  cerr << "Error on writing buffer to " << options.anova_indexfile << '\n';
	  exit(errno);
	}
    }
  else
    {
      FILE * a_fh = fopen(options.anovafile.c_str(),"a"),*id_fh=NULL;
      int a_fd = fileno(a_fh),id_fd=-1;
      ostringstream buffer;
      buffer << options.record_no << ' ' << lseek( a_fd, 0, SEEK_CUR );
      if (!options.idfile.empty())
	{
	  id_fh = fopen(options.idfile.c_str(),"a");
	  id_fd = fileno(id_fh);
	  buffer <<' ' << lseek(id_fd,0,SEEK_CUR) << '\n';
	  if(::write(id_fd,
		     idbuffer.str().c_str(),
		     idbuffer.str().size()) == -1 )
	    {
	      cerr << "Error on writing buffer to " << options.idfile << '\n';
	      exit(10);
	    }
	  fclose(id_fh);
	}
      else
	{
	  buffer << '\n';
	}
      if( ::write(ai_fd,buffer.str().c_str(),buffer.str().size()) == -1 )
	{
	  cerr << "Error on writing buffer to " << options.anova_indexfile << '\n';
	  exit(errno);
	}
      buffer.str(string());
      
      //write out the CC buffer
      if(::write( a_fd, ccbuffer.str().c_str(), ccbuffer.str().size() ) == -1 )
	{
	  cerr << "Error on writing buffer to " << options.anovafile << '\n';
	  exit(errno);
	}
      fflush(a_fh);
      fclose(a_fh);
    }

  // FILE * a_fh = fopen(options.anovafile.c_str(),"a");
  // int a_fd = fileno(a_fh);
  // struct flock a_lock = get_whole_flock();
  // assert( a_lock.l_type == F_WRLCK );
  // assert( a_lock.l_whence == SEEK_SET );
  // assert( a_lock.l_start == 0 );
  // assert( a_lock.l_len == 0 );
  // if (fcntl(a_fd,F_SETLKW,&a_lock) == -1)
  //   {
  //     cerr << "ERROR: could not obtain lock on " << options.anovafile << '\n';
  //     exit(10);
  //   }
  
  // ostringstream buffer;
  // buffer << options.record_no << ' ' << lseek( a_fd, 0, SEEK_CUR ) << '\n';
  // if( ::write(ai_fd,buffer.str().c_str(),buffer.str().size()) == -1 )
  //   {
  //     cerr << "Error on writing buffer to " << options.anova_indexfile << '\n';
  //     exit(errno);
  //   }
  // buffer.str(string());
  
  // //write out the CC buffer
  // if(::write( a_fd, ccbuffer.str().c_str(), ccbuffer.str().size() ) == -1 )
  //   {
  //     cerr << "Error on writing buffer to " << options.anovafile << '\n';
  //     exit(errno);
  //   }
  
  // //Unlock files
  // a_lock.l_type = F_UNLCK;
  // if (fcntl(a_fd,F_UNLCK,&a_lock) == -1)
  //   {
  //     cerr << "ERROR: could not relesae lock on " << options.anovafile << '\n';
  //     exit(10);
  //   }
  // fflush(a_fh);
  // fclose(a_fh);

  //release lock on index file
  /*
  ai_lock.l_type = F_UNLCK;
  if (fcntl(ai_fd,F_UNLCK,&ai_lock) == -1)
    {
      cerr << "ERROR: could not release lock on " << options.anova_indexfile << '\n';
      exit(10);
    }
  */
  fflush(ai_fh);
  fclose(ai_fh);
  ai_lock.unlock();
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
    ("ccfile,c",value<string>(&rv.anovafile)->default_value(string()),"File to write case/control data")
    ("ccindex,I",value<string>(&rv.anova_indexfile)->default_value(string()),"File to write index info for case/control data file")
    ("recordno,r",value<unsigned>(&rv.record_no)->default_value(0),"Record number to look up in indexfile")
    ("maxcases,n",value<unsigned>(&rv.ncases)->default_value(0),"Maximum number of cases to sample")
    ("maxcontrols,N",value<unsigned>(&rv.ncontrols)->default_value(0),"Maximum Number of controls to sample")
    ("seed,S",value<unsigned>(&rv.seed)->default_value(0),"Random number seed")
    ("threshold,t",value<double>(&rv.case_proportion)->default_value(-1.),"Proportion of population to be labelled as putative cases.  Value must be 0 < t < 1. E.g., individuals with phenotypic values larger than the (1-t)th quantile of phenotypic values in the entire population are potential cases.  For example, in Thornton, Foran, and Long (2013), we used t = 0.15, meaning that the upper 15% of phenotypic values were treated as possible cases.")
    ("control-range",value<double>(&rv.crange)->default_value(0.5),"A putatitve control is defined as the population mean phenotype +/- control-range*SD, where SD is the standard deviation of the population phenotype.  Default is what was used in Thornton, Foran, and Long (2013)")
    ("gzin","Input is gzipped.  Default is to assume uncompressed input.")
    ("gzout","Write output as gzipped file.  Default is to write uncompressed output.")
    ("ids",value<string>(&rv.idfile)->default_value(string()),"Write individual identifiers to output file name.  Useful of you want to go back to the raw genotype info in the main population file later on.")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if(argc == 1 || vm.count("help"))
    {
      cerr << desc << '\n';
      exit(0);
    }

  if (rv.files_undef())
    {
      rv.report_empty(cerr);
      exit(10);
    }

  if (! rv.params_ok() )
    {
      exit(10);
    }

  if( vm.count("gzin") )
    {
      rv.gzinput=true;
    }
  if(vm.count("gzout"))
    {
      rv.gzoutput = true;
    }
  return rv;
}

bool params::files_undef( void ) const
{
  return ( indexfile.empty() ||
	   popfile.empty() ||
	   phenofile.empty() ||
	   anovafile.empty() ||
	   anova_indexfile.empty() );
}

void params::report_empty( std::ostream & out ) const
{
  if( indexfile.empty() )
    {
      out << "Error: index file for input data not defined.  Use -i option.\n";
    }
  if( popfile.empty() )
    {
      out << "Error: Population data file for input data not defined.  Use -p option.\n";
    }
  if( phenofile.empty() )
    {
      out << "Error: Phenotype data file for input data not defined.  Use -P option.\n";
    }
  if( anovafile.empty() )
    {
      out << "Error: Output data file for case/control genotypes not specified.  Use -c option.\n";
    }
  if( anova_indexfile.empty() )
    {
      out << "Error: Output file name for index file not specified.  Use -I option.\n";
    }
}

bool params::params_ok( void ) const
{
  bool ok = true;

  if ( ! ncases )
    {
      ok = false;
      cerr << "Error: # of cases = " << ncases << '\n';
    }
  if ( ! ncontrols )
    {
      ok = false;
      cerr << "Error: # of controls = " << ncontrols << '\n';
    }

  if ( ! (case_proportion > 0. && case_proportion < 1.) )
    {
      ok = false;
      cerr << "Error: proportion of cases must be 0 < x < 1.  Input value was " << case_proportion << '\n';
    }

  if ( crange <= 0. )
    {
      ok = false;
      cerr << "Error: bounds on controls (--control-range) must be positive (> 0.).  You entered " << crange << '\n';
    }

  ifstream in(indexfile.c_str());
  if(!in)
    {
      cerr << "Error, could not open " << indexfile << " for reading\n";
      ok = false;
    }
  in.close();

  in.open(popfile.c_str());
  if(!in)
    {
      cerr << "Error, could not open " << popfile << " for reading\n";
      ok = false;
    }
  in.close();

  in.open(phenofile.c_str());
  if(!in)
    {
      cerr << "Error, could not open " << phenofile << " for reading\n";
      ok = false;
    }
  in.close();

  if (! ok )
    {
      cerr << "Please use -h option to see help\n";
    }
  return ok;
}

std::pair<double,double> phenosums(const vector<pair<double,double> > & phenos, const double & case_proportion, double * cutoff)
{
  vector<double> pcopy;
  for(unsigned i=0;i<phenos.size();++i)
    {
      pcopy.push_back(phenos[i].first+phenos[i].second);
    }
  sort(pcopy.begin(),pcopy.end());
  //get the upper case_proportion-th'd quantile from the pheno dist
  *cutoff = gsl_stats_quantile_from_sorted_data(&pcopy[0],1,pcopy.size(),1.-case_proportion);
  return std::make_pair( gsl_stats_mean(&pcopy[0],1,pcopy.size()),
			 gsl_stats_sd(&pcopy[0],1,pcopy.size()) );
}
