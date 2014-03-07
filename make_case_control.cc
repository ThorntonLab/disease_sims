/*
  Make case/control "panel" assuming a liability threshold on phenotype.
  In this case, individuals with phenotypes >= the x-th quantile
  of the population's phenotypic distribution are potential cases.

  There are two output files:

  1. anovafile_index:  this file records the record_no passed to the command line,
  and the offset (bytes) of where that record begins in anovafile

  2. anovafile: This fine is written in binary format.  It contains a 
  series of records of case/control panels.  Each record has the following format:
  Four integers: number of cases, number of controls, NN = number of neutral mutations, NC = number of causative mutations
  NN + NC doubles, which are the positions of the mutations (neutral then causative mutations, respectively).
  Then, for each case, and then each control, there are NN + NC short signed integers
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
//#include <locking_routines.hpp>

#include <Sequence/PolyTableFunctions.hpp>

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

using namespace std;
using namespace Sequence;
using namespace KTfwd;

typedef gamete_base<TFLmtype> gtype;
typedef list<TFLmtype> mlist;
typedef list<gtype> glist;

struct params
{
  string indexfile,popfile,phenofile,anovafile,anovafile_index;
  unsigned twoN,record_no,ncases,ncontrols,seed;
  double case_proportion;
};

params process_command_line(int argc, char ** argv);


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
  vector< double > phenotypes;
  unsigned nphenos;
  phenostream.read( reinterpret_cast<char *>(&nphenos), sizeof(unsigned) );
  for( unsigned i = 0 ; i < nphenos ; ++i )
    {
      //x is the genetic contribution to phenotype. y is the Gaussian noise from the simulation.
      //Phenotype of the individual is x+y
      double x,y;
      phenostream.read( reinterpret_cast<char *>(&x), sizeof(double) );
      phenostream.read( reinterpret_cast<char *>(&y), sizeof(double) );
      phenotypes.push_back( x+y );
    }
  phenostream.close();
  //obtain file lock on index ASAP
  // FILE * ai_fh = fopen(anova_indexfile,"a");
  // int ai_fd = fileno(ai_fh);

  // flock ai_lock = get_whole_flock();
  // //make sure our locking functions work...
  // assert( ai_lock.l_type == F_WRLCK );
  // assert( ai_lock.l_whence == SEEK_SET );
  // assert( ai_lock.l_start == 0 );
  // assert( ai_lock.l_len == 0 );
  // if (fcntl(ai_fd,F_SETLKW,&ai_lock) == -1)
  //   {
  //     cerr << "ERROR: could not obtain lock on " << anova_indexfile << '\n';
  //     exit(10);
  //   }

  // FILE * a_fh = fopen(anovafile,"a");
  // int a_fd = fileno(a_fh);
  // flock a_lock = get_whole_flock();
  // assert( a_lock.l_type == F_WRLCK );
  // assert( a_lock.l_whence == SEEK_SET );
  // assert( a_lock.l_start == 0 );
  // assert( a_lock.l_len == 0 );
  // if (fcntl(a_fd,F_SETLKW,&a_lock) == -1)
  //   {
  //     cerr << "ERROR: could not obtain lock on " << anovafile << '\n';
  //     exit(10);
  //   }

  // //buffer and write output

  // //first, the index info.
  // ostringstream buffer;
  // buffer << record_no << ' ' << lseek( a_fd, 0, SEEK_CUR ) << '\n';
  // write(ai_fd,buffer.str().c_str(),buffer.str().size());
  // buffer.str(string());

  // //now, the data that will be processed later for actual GWAS
  // buffer.write(reinterpret_cast<const char *>(&ncontrols),sizeof(unsigned));
  // buffer.write(reinterpret_cast<const char *>(&ncases),sizeof(unsigned));
  // unsigned neutral_muts = neutral.numsites(),causative_muts=causative.numsites();
  // buffer.write(reinterpret_cast<char *>(&neutral_muts),sizeof(unsigned));
  // buffer.write(reinterpret_cast<char *>(&causative_muts),sizeof(unsigned));

  // for(SimData::const_pos_iterator i = neutral.pbegin() ; i != neutral.pend() ; ++i)
  //   {
  //     double pos = *i;
  //     buffer.write(reinterpret_cast<char *>(&pos),sizeof(double));
  //   }
  // for(SimData::const_pos_iterator i = causative.pbegin() ; i != causative.pend() ; ++i)
  //   {
  //     double pos = *i;
  //     buffer.write(reinterpret_cast<char *>(&pos),sizeof(double));
  //   }
  // write(a_fd,buffer.str().c_str(),buffer.str().size());
  // buffer.str(string());

  // vector<char> minors_neutral = get_minors(neutral);
  // assert(minors_neutral.size() == neutral.numsites());
  // vector<char> minors_causative = get_minors(causative);
  // assert(minors_causative.size() == causative.numsites());

  // //write the data in binary to output file
  // write_data(a_fd, 0, 2*(ncontrols+ncases),
  // 	     neutral,causative,
  // 	     minors_neutral,minors_causative,
  // 	     phenodata);

  // a_lock.l_type = F_UNLCK;
  // if (fcntl(a_fd,F_UNLCK,&a_lock) == -1)
  //   {
  //     cerr << "ERROR: could not obtain lock on " << anovafile << '\n';
  //     exit(10);
  //   }
  // fclose(a_fh);
  // ai_lock.l_type = F_UNLCK;
  // if (fcntl(ai_fd,F_UNLCK,&ai_lock) == -1)
  //   {
  //     cerr << "ERROR: could not release lock on " << anova_indexfile << '\n';
  //     exit(10);
  //   }
  // fclose(ai_fh);
}

// void get_pheno_stats(const population & p,const double & case_proportion,
// 		     double * pheno_cutoff,double * mean_pheno,
// 		     double * sd_pheno)
// {
//   vector<double> phenos;
//   for(unsigned i = 0 ; i < p.phenotypes.size() ; ++i )
//     {
//       phenos.push_back(p.phenotypes[i].first + p.phenotypes[i].second);
//     }
//   sort(phenos.begin(),phenos.end());
//   *pheno_cutoff = gsl_stats_quantile_from_sorted_data(&phenos[0],1,
// 						      phenos.size(),
// 						      1. - case_proportion);
//   *mean_pheno = gsl_stats_mean(&phenos[0],1,phenos.size());
//   *sd_pheno = gsl_stats_sd(&phenos[0],1,phenos.size());
// }


params process_command_line(int argc, char ** argv)
{
  params rv;

  return rv;
}
