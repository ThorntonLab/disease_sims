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
  Then, for each control, and then each case, there are NN + NC short signed integers
  corresponding to the genotypes and neutral and causative sites, resp. 

  A genotype can take on 1 of 3 values, and is equal to the # of copies of the minor allele:
  0 = homozygote for major allele
  1 = heterozygote
  2 = homozygote for minor allele

  The genotypes are stored in the same order as the positions.

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
#include <locking_routines.hpp>

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

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <fwdpp/IO.hpp>
#include <fwdpp/sampling_functions.hpp>


using namespace std;
using namespace Sequence;
using namespace KTfwd;

struct params
{
  string indexfile,popfile,phenofile,anovafile,anova_indexfile;
  unsigned twoN,record_no,ncases,ncontrols,seed;
  double case_proportion;
};

params process_command_line(int argc, char ** argv);

std::pair<double,double> phenosums(const vector<pair<double,double> > & phenos, const double & case_proportion, double * cutoff);

int main(int argc, char ** argv)
{
  params options = process_command_line(argc,argv);
  int argn = 1;
  // const char * indexfile = argv[argn++];
  // const char * hapfile = argv[argn++];
  // const char * phenofile = argv[argn++];
  // const unsigned & twoN = atoi(argv[argn++]);
  // const unsigned & record_no = atoi(argv[argn++]);
  // const char * anovafile = argv[argn++];
  // const char * anova_indexfile = argv[argn++];
  // const double case_proportion = atof(argv[argn++]); //e.g., 0.1 means the top 10% of the population's phenotypes are considered "diseased"
  // const unsigned ncases = atoi(argv[argn++]);
  // const unsigned ncontrols = atoi(argv[argn++]);
  // const unsigned seed = atoi(argv[argn++]);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,options.seed);
  
  //Read in the population
  glist gametes;
  mlist mutations;
  vector< pair< glist::iterator,glist::iterator > > diploids;

  long pop_offset,pheno_offset,effect_offset;
  unsigned ith_rep;
  bool found = true;
  //get offset of the population and the phenotypes from the indexfile
  ifstream index( options.indexfile.c_str() );
  while( !found && !index.eof() )
    {
      index >> ith_rep >> effect_offset >> pheno_offset >> pop_offset >> ws;
      if( ith_rep == options.record_no )
	{
	  found = true;
	}
    }

  if ( ! found )
    {
      cerr << "Error: replicate number " << options.record_no << " not found in " << options.indexfile << '\n';
      exit(10);
    }
  index.close();

  ifstream popstream( options.popfile.c_str() );
  popstream.seekg( pop_offset );
  read_binary_pop( &gametes, &mutations, &diploids, boost::bind(mreader(),_1),popstream );
  popstream.close();

  //Read in the phenotype data
  ifstream phenostream( options.phenofile.c_str() );
  phenostream.seekg( pheno_offset );
  vector< pair<double,double> > phenotypes;
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
      else if ( P >= mean_sd.first - mean_sd.second &&
		P <= mean_sd.first + mean_sd.second )
	{
	  put_controls.push_back(i);
	}
    }

  //Randomize lists just for fun
  boost::function< size_t (size_t) > rand = boost::bind(&gsl_ran_flat, r, 0,double(put_cases.size()));
  random_shuffle(put_cases.begin(),put_cases.end(),rand);
  rand = boost::bind(&gsl_ran_flat, r, 0,double(put_controls.size()));
  random_shuffle(put_controls.begin(),put_controls.end(),rand);

  cc_intermediate * ccblocks = new cc_intermediate;

  *ccblocks = process_population(diploids,phenotypes,
				 put_controls,
				 put_cases,
				 options.ncontrols,
				 options.ncases);
  //free up memory
  diploids.clear();
  gametes.clear();
  mutations.clear();

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
      //neutral genotypes for this individual
      for( unsigned site = 0 ; site < ccblocks->neutral.numsites() ; ++site )
	{
	  //count # copies of minor allele at this site in this individual
	  unsigned cminor = ( (ccblocks->neutral[ind][site] == ccblocks->min_n[site]) ? 1 : 0 ) +
	    ( (ccblocks->neutral[ind+1][site] == ccblocks->min_n[site]) ? 1 : 0 );
	  ccbuffer.write( reinterpret_cast<char *>(&cminor), sizeof(unsigned) );
	}
      //causative genotypes for this individual
      for( unsigned site = 0 ; site < ccblocks->causative.numsites() ; ++site )
	{
	  //count # copies of minor allele at this site in this individual
	  unsigned cminor = ( (ccblocks->causative[ind][site] == ccblocks->min_n[site]) ? 1 : 0 ) +
	    ( (ccblocks->causative[ind+1][site] == ccblocks->min_n[site]) ? 1 : 0 );
	  ccbuffer.write( reinterpret_cast<char *>(&cminor), sizeof(unsigned) );
	}
    }

  //Now, output # of causative mutations on each haplotype carried by this diploid
  for( unsigned ind = 0 ; ind < ccblocks->neutral.size() ; ind += 2 ) 
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

   flock ai_lock = get_whole_flock();

   //make sure our locking functions work...
   assert( ai_lock.l_type == F_WRLCK );
   assert( ai_lock.l_whence == SEEK_SET );
   assert( ai_lock.l_start == 0 );
   assert( ai_lock.l_len == 0 );
   if (fcntl(ai_fd,F_SETLKW,&ai_lock) == -1)
     {
       cerr << "ERROR: could not obtain lock on " << options.anova_indexfile << '\n';
       exit(10);
     }

   FILE * a_fh = fopen(options.anovafile.c_str(),"a");
   int a_fd = fileno(a_fh);
   flock a_lock = get_whole_flock();
   assert( a_lock.l_type == F_WRLCK );
   assert( a_lock.l_whence == SEEK_SET );
   assert( a_lock.l_start == 0 );
   assert( a_lock.l_len == 0 );
   if (fcntl(a_fd,F_SETLKW,&a_lock) == -1)
     {
       cerr << "ERROR: could not obtain lock on " << options.anovafile << '\n';
       exit(10);
     }

   ostringstream buffer;
   buffer << options.record_no << ' ' << lseek( a_fd, 0, SEEK_CUR ) << '\n';
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
   
   //Unlock files
   a_lock.l_type = F_UNLCK;
   if (fcntl(a_fd,F_UNLCK,&a_lock) == -1)
     {
       cerr << "ERROR: could not relesae lock on " << options.anovafile << '\n';
       exit(10);
     }
   fflush(a_fh);
   fclose(a_fh);
   ai_lock.l_type = F_UNLCK;
   if (fcntl(ai_fd,F_UNLCK,&ai_lock) == -1)
     {
       cerr << "ERROR: could not release lock on " << options.anova_indexfile << '\n';
       exit(10);
      }
   fflush(ai_fh);
   fclose(ai_fh);
}

params process_command_line(int argc, char ** argv)
{
  params rv;

  return rv;
}

std::pair<double,double> phenosums(const vector<pair<double,double> > & phenos, const double & case_proportion, double * cutoff)
{
  vector<double> pcopy;
  for(unsigned i=0;i<phenos.size();++i)
    {
      pcopy.push_back(phenos[i].first+phenos[i].second);
    }
  sort(pcopy.begin(),pcopy.end());
  return std::make_pair( gsl_stats_mean(&pcopy[0],1,pcopy.size()),
			 gsl_stats_sd(&pcopy[0],1,pcopy.size()) );
  //get the upper case_proportion-th'd quantile from the pheno dist
  *cutoff =  gsl_stats_quantile_from_sorted_data(&pcopy[0],1,pcopy.size(),1.-case_proportion);
}
