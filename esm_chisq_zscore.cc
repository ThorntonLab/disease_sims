/*
  Outputs the ESM statistic from Thornton, Foran, and Long (2013) PLoS Genetics.

  Instead of Fisher's Exact Test (FET) as done in Thornton et al., we use the 
  chi-squared statistic here, with a continuity correction applied.  The chi-squared
  is much, much faster to compute.

  The output of the program is two numbers.  The first is the ID number of the
  data set analyzed.  The second number is a z-score obtained based on the observed ESM
  statistic and P permutations of the case/control data.
*/

#include <readCC.hpp>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <numeric>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include <gwas_stats.hpp>
#include <locking_routines.hpp>

using namespace std;

void permute( const CCblock * ccdata,
	      const unsigned & K,
	      const vector<short> * keep,
	      const unsigned & nperms,
	      const unsigned & seed,
	      vector<double> * rv );

int main(int argc, char ** argv)
{
  int argn = 1;
  const char * ccindexfile = argv[argn++];
  const char * ccfile = argv[argn++];
  const unsigned record_no = atoi(argv[argn++]);
  const double minfreq = atof(argv[argn++]);
  const double maxfreq = atof(argv[argn++]);
  const unsigned K = atoi(argv[argn++]);
  const double rsq_cutoff = atof(argv[argn++]);
  const unsigned max_controls = atoi(argv[argn++]);
  const unsigned max_cases = atoi(argv[argn++]);
  const unsigned nperms = atoi(argv[argn++]);
  const char * outfile = argv[argn++];
  const unsigned nthreads = atoi(argv[argn++]);
  const unsigned seed = atoi(argv[argn++]);

  bool fail = true;
  CCblock ccdata = read_CC_record(ccindexfile,ccfile,record_no,&fail,max_controls,max_cases,true);
  if(fail)
    {
      cerr << "Error: could not find record " << record_no << " in " << ccindexfile << '\n';
      exit(10);
    }

  vector<short> ccstatus(ccdata.ncontrols,0);
  for(unsigned i=0;i<ccdata.ncases;++i)
    {
      ccstatus.push_back(1);
    }
  assert(ccstatus.size() == ccdata.geno_matrix[0].size());
  vector<short> keep = filter_sites(ccdata,minfreq,maxfreq,rsq_cutoff);

  vector<double> chisqs = chisq_per_marker( &ccdata, &keep, ccstatus );

  double obs_esm = esm(chisqs,K);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);

  vector<double> permstats;
  if(nthreads == 1)
    {
      permute(&ccdata,K,&keep,nperms,gsl_rng_get(r),&permstats);
    }
  else if (nthreads > 1)
    {
      vector<vector<double> > permstats_thread(nthreads,vector<double>());
      unsigned nperms_thread = nperms/nthreads;
      boost::thread_group tg;
      for( unsigned i=0;i<nthreads;++i)
	{
	  if( i == 0 && nperms_thread*nthreads < nperms )
	    {
	      tg.add_thread( new boost::thread( boost::bind(permute,&ccdata,K,&keep,nperms_thread + (nperms-nperms_thread*nthreads),gsl_rng_get(r),&permstats_thread[i]) ) );
	    }
	  else
	    {
	      tg.add_thread( new boost::thread( boost::bind(permute,&ccdata,K,&keep,nperms_thread,gsl_rng_get(r),&permstats_thread[i]) ) );
	    }
	}
      tg.join_all();
      
      
      for(unsigned i=0;i<permstats_thread.size();++i)
	{
	  assert( permstats_thread[i].size() >= nperms_thread );
	  copy(permstats_thread[i].begin(),permstats_thread[i].end(),
	       std::back_inserter(permstats));
	}
      assert( permstats.size() == nperms );
      permstats_thread.clear();

    }
  sort(permstats.begin(),permstats.end());
  double perm_p = double( count_if(permstats.begin(),permstats.end(),
				   boost::bind(greater<double>(),_1,obs_esm ) ) )/double(nperms);
  double mean = gsl_stats_mean(&permstats[0],1,permstats.size());
  double sd = gsl_stats_sd(&permstats[0],1,permstats.size());
  double z = (obs_esm-mean)/sd;

  FILE * fp = fopen(outfile,"a");
  int fd = fileno(fp);

  flock fd_lock = get_whole_flock();

  if (fcntl(fd,F_SETLKW,&fd_lock) == -1)
    {
      cerr << "ERROR: could not obtain lock on " << outfile << '\n';
      exit(10);
    }

  fprintf(fp,"%u\t%lf\t%lf\n",record_no,perm_p,z);

  fd_lock.l_type = F_UNLCK;
  if (fcntl(fd,F_UNLCK,&fd_lock) == -1)
    {
      cerr << "ERROR: could not release lock on " << outfile << '\n';
      exit(10);
    }
  fclose(fp);
}

void permute( const CCblock * ccdata,
	      const unsigned & K,
	      const vector<short> * keep,
	      const unsigned & nperms,
	      const unsigned & seed,
	      vector<double> * rv )
{
  vector<short> ccstatus(ccdata->ncontrols,0);
  for(unsigned i=0;i<ccdata->ncases;++i)
    {
      ccstatus.push_back(1);
    }
  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);
  boost::function< size_t (size_t) > rand = boost::bind(&gsl_ran_flat, r, 0,double(ccstatus.size()));
  for(unsigned i=0;i<nperms;++i)
    {
      random_shuffle( ccstatus.begin(), ccstatus.end(), rand );
      vector<double> chisqs = chisq_per_marker( ccdata, keep, ccstatus );
      rv->push_back( esm(chisqs,K) );
    }
  gsl_rng_free(r);
}
