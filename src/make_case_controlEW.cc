/*
  Make case/control "panel" assuming a liability threshold on phenotype.
  In this case, individuals with phenotypes >= the x-th quantile
  of the population's phenotypic distribution are potential cases.

  This program differs from make_case_control in that the input
  are assumed to come from a simulation done under the "EW2010" model.

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

#include <diseaseSims/mutation_with_age.hpp>
#include <diseaseSims/ccintermediate.hpp>
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
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/program_options.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <fwdpp/IO.hpp>
#include <fwdpp/sampling_functions.hpp>

#include <readSimOutput.hpp>

using namespace std;
using namespace Sequence;
using namespace KTfwd;
using namespace boost::interprocess;
using namespace boost::program_options;

struct params
{
  string indexfile,popfile,anovafile,anova_indexfile,idfile;
  unsigned N,record_no,ncases,ncontrols,seed;
  double case_proportion,crange,tau,sigma,pdelta;
  params();
  bool files_undef( void ) const;
  void report_empty( std::ostream & out ) const;
  bool params_ok( void ) const;
};

params::params() : indexfile(string()),
		   popfile(string()),
		   anovafile(string()),
		   anova_indexfile(string()),
		   idfile(string()),
		   N(0),
		   record_no(std::numeric_limits<unsigned>::max()),
		   ncases(0),
		   ncontrols(0),
		   seed(0),
		   case_proportion(1.),
		   crange(0.5),
		   tau(1.),
		   sigma(1.),
		   pdelta(0.5)
{
}

params process_command_line(int argc, char ** argv);

vector<pair<double,double> > EWphenos( const vector< pair< glist::iterator,glist::iterator > > & diploids,
				       const vector< pair<mlist::iterator, double> > & ESIZES );

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

  gzFile gzin = gzopen(options.popfile.c_str(),"rb");
  gzseek( gzin, index.hoffset(options.record_no), 0);
  read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
  gzclose(gzin);

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

  //Now, we need to assign effect sizes to each mutation that is segregating.
  vector< pair<mlist::iterator, double> > ESIZES;
  for( auto mitr = mutations.begin(); mitr != mutations.end() ; ++mitr )
    {
      if( ! mitr->neutral )
	{
	  double __delta = (gsl_rng_uniform(r) <= options.pdelta) ? 1. : -1.;
	  ESIZES.push_back( make_pair(mitr,__delta*pow(4.*double(options.N)*(mitr->s),options.tau)*(1. + gsl_ran_gaussian(r,options.sigma))) );
	}
    }

  //Generate the phenotype data
  vector< pair<double,double> > phenotypes = EWphenos(diploids,ESIZES);
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
  vector< unsigned > put_cases,put_controls;
  grab_putative_CC(mean_sd,phenotypes,options.crange,cutoff,put_controls,put_cases);


  /*
    Check: do putative cases and controls overlap?

    If so, can we get rid of the overlap and still have sufficient #s in each category?

    If not, we warn that the CC data may be invalid, as cases and controls will share individuals.
  */

  vector< unsigned > isect;
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
      vector< unsigned > ucontrols, ucases;
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
  random_shuffle( put_cases.begin(),put_cases.end(),[&r](int n) { return gsl_rng_get(r) % n; });
  random_shuffle( put_controls.begin(),put_controls.end(),[&r](int n) { return gsl_rng_get(r) % n; });

#ifndef NDEBUG
  for ( auto x : put_cases ) {
    assert( x < diploids.size() );
  }
  for ( auto x : put_controls ) {
    assert( x < diploids.size() );
  }
#endif
  cc_intermediate ccblocks(process_population(diploids,phenotypes,
					      put_controls,
					      put_cases,
					      options.ncontrols,
					      options.ncases) );
  assert( ccblocks.min_n.size() == ccblocks.neutral.numsites() );
  assert( ccblocks.min_c.size() == ccblocks.causative.numsites() );
  
  //3. Output to buffers

  //first, the ccdata to a buffer
  ostringstream ccbuffer;
  ccbuffer << ccblocks;

  ostringstream idbuffer;
  if( !options.idfile.empty() ) //then we want to write the individual id indexes
    {
      idbuffer.write( reinterpret_cast<char *>(&ccblocks.ncontrols), sizeof(unsigned) );  
      idbuffer.write( reinterpret_cast<char *>(&ccblocks.ncases), sizeof(unsigned) );  
      idbuffer.write( reinterpret_cast<char *>(&ccblocks.control_ids[0]), ccblocks.ncontrols*sizeof(unsigned));
      idbuffer.write( reinterpret_cast<char *>(&ccblocks.case_ids[0]), ccblocks.ncases*sizeof(unsigned));
    }

  //4. Output to files

  //obtain file lock on index ASAP
  FILE * ai_fh = fopen(options.anova_indexfile.c_str(),"a");
  int ai_fd = fileno(ai_fh);
  
  file_lock ai_lock(options.anova_indexfile.c_str());
  scoped_lock<file_lock> s_lock(ai_lock);

  gzFile gzout = gzopen( options.anovafile.c_str(),"a" );
  int written = gzwrite( gzout, ccbuffer.str().c_str(), ccbuffer.str().size() );
  if(!written)
    {
      cerr << "Error writing to " << options.anovafile << '\n';
      exit(10);
    }
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

  fflush(ai_fh);
  ai_lock.unlock();
  fclose(ai_fh);
  exit(0);
}

/*
  This is likely not the fastest way to do things...

  ...screw it.
*/
vector<pair<double,double> > EWphenos( const vector< pair< glist::iterator,glist::iterator > > & diploids,
				       const vector< pair<mlist::iterator, double> > & ESIZES )
{
  vector<pair<double,double> > rv;

  for( unsigned i = 0 ; i < diploids.size() ; ++i )
    {
      double trait = 0.;
      for( const auto & mitr : diploids[i].first->smutations )
	{
	  auto __cc = find_if( ESIZES.begin(), ESIZES.end(),
			       [&mitr]( const pair<mlist::iterator, double> & __p ) 
			       {
				 return __p.first == mitr;
			       } );
	  assert(__cc != ESIZES.end());
	  if(__cc == ESIZES.end())
	    {
	      cerr << "FATAL ERROR: mutation carried by a gamete is not present in the population. "
		   << "Line " << __LINE__ << " of " << __FILE__ << '\n';
	      exit(EXIT_FAILURE);
	    }
	  trait += __cc->second;
	}
      for( const auto & mitr : diploids[i].second->smutations )
	{
	  auto __cc = find_if( ESIZES.begin(), ESIZES.end(),
			       [&mitr]( const pair<mlist::iterator, double> & __p ) 
			       {
				 return __p.first == mitr;
			       } );
	  assert(__cc != ESIZES.end());
	  if(__cc == ESIZES.end())
	    {
	      cerr << "FATAL ERROR: mutation carried by a gamete is not present in the population. "
		   << "Line " << __LINE__ << " of " << __FILE__ << '\n';
	      exit(EXIT_FAILURE);
	    }
	  trait += __cc->second;
	}
      rv.emplace_back( make_pair( trait, 0. ) ) ; //It's all "G", kids!
    }
  return rv;
}

params process_command_line(int argc, char ** argv)
{
  params rv;

  options_description desc("Process output from TFL2013_ind and generate case/control panel");
  desc.add_options()
    ("help,h", "Produce help message")
    ("indexfile,i",value<string>(&rv.indexfile)->default_value(string()),"Index file output by TFL2013_ind")
    ("popfile,p",value<string>(&rv.popfile)->default_value(string()),"Population file output by TFL2013_ind")
    ("popsize",value<unsigned>(&rv.N)->default_value(0u),"The population size used in the simulation to convert S to s")
    ("tau,T",value<double>(&rv.tau)->default_value(1.),"The power applied to S")
    ("sigma",value<double>(&rv.sigma)->default_value(1.),"The std. deviation use on e")
    ("symmetry",value<double>(&rv.pdelta)->default_value(0.5),"Probability that delta = 1, otherwise delta = -1")
    ("ccfile,c",value<string>(&rv.anovafile)->default_value(string()),"File to write case/control data")
    ("ccindex,I",value<string>(&rv.anova_indexfile)->default_value(string()),"File to write index info for case/control data file")
    ("recordno,r",value<unsigned>(&rv.record_no)->default_value(0),"Record number to look up in indexfile")
    ("maxcases,n",value<unsigned>(&rv.ncases)->default_value(0),"Maximum number of cases to sample")
    ("maxcontrols,N",value<unsigned>(&rv.ncontrols)->default_value(0),"Maximum Number of controls to sample")
    ("seed,S",value<unsigned>(&rv.seed)->default_value(0),"Random number seed")
    ("threshold,t",value<double>(&rv.case_proportion)->default_value(-1.),"Proportion of population to be labelled as putative cases.  Value must be 0 < t < 1. E.g., individuals with phenotypic values larger than the (1-t)th quantile of phenotypic values in the entire population are potential cases.  For example, in Thornton, Foran, and Long (2013), we used t = 0.15, meaning that the upper 15% of phenotypic values were treated as possible cases.")
    ("control-range",value<double>(&rv.crange)->default_value(0.5),"A putatitve control is defined as the population mean phenotype +/- control-range*SD, where SD is the standard deviation of the population phenotype.  Default is what was used in Thornton, Foran, and Long (2013)")
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

  return rv;
}

bool params::files_undef( void ) const
{
  return ( indexfile.empty() ||
	   popfile.empty() ||
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

  if (! ok )
    {
      cerr << "Please use -h option to see help\n";
    }
  return ok;
}


