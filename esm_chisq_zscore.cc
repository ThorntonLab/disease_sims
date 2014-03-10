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

#include <boost/program_options.hpp>
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

struct params
{
  string ccindexfile,ccfile,outfile;
  unsigned record_no,K,max_controls,max_cases,nperms,nthreads,seed;
  double minfreq,maxfreq,rsq_cutoff;
};

params parse_command_line(int argc, char ** argv);

int main(int argc, char ** argv)
{
  // int argn = 1;
  // //const char * ccindexfile = argv[argn++];
  // //const char * ccfile = argv[argn++];
  // const unsigned record_no = atoi(argv[argn++]);
  // const double minfreq = atof(argv[argn++]);
  // const double maxfreq = atof(argv[argn++]);
  // const unsigned K = atoi(argv[argn++]);
  // const double rsq_cutoff = atof(argv[argn++]);
  // const unsigned max_controls = atoi(argv[argn++]);
  // const unsigned max_cases = atoi(argv[argn++]);
  // const unsigned nperms = atoi(argv[argn++]);
  // //const char * outfile = argv[argn++];
  // const unsigned nthreads = atoi(argv[argn++]);
  // const unsigned seed = atoi(argv[argn++]);

  params options = parse_command_line(argc, argv);

  bool fail = true;
  CCblock ccdata = read_CC_record(options.ccindexfile.c_str(),options.ccfile.c_str(),options.record_no,&fail,options.max_controls,options.max_cases,true);
  if(fail)
    {
      cerr << "Error: could not find record " << options.record_no << " in " << options.ccindexfile << '\n';
      exit(10);
    }

  vector<short> ccstatus(ccdata.ncontrols,0);
  for(unsigned i=0;i<ccdata.ncases;++i)
    {
      ccstatus.push_back(1);
    }
  assert(ccstatus.size() == ccdata.geno_matrix[0].size());
  vector<short> keep = filter_sites(ccdata,options.minfreq,options.maxfreq,options.rsq_cutoff);

  vector<double> chisqs = chisq_per_marker( &ccdata, &keep, ccstatus );

  double obs_esm = esm(chisqs,options.K);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,options.seed);

  vector<double> permstats;
  if(options.nthreads == 1)
    {
      permute(&ccdata,options.K,&keep,options.nperms,gsl_rng_get(r),&permstats);
    }
  else if (options.nthreads > 1)
    {
      vector<vector<double> > permstats_thread(options.nthreads,vector<double>());
      unsigned nperms_thread = options.nperms/options.nthreads;
      boost::thread_group tg;
      for( unsigned i=0;i<options.nthreads;++i)
	{
	  if( i == 0 && nperms_thread*options.nthreads < options.nperms )
	    {
	      tg.add_thread( new boost::thread( boost::bind(permute,&ccdata,options.K,&keep,nperms_thread + (options.nperms-nperms_thread*options.nthreads),gsl_rng_get(r),&permstats_thread[i]) ) );
	    }
	  else
	    {
	      tg.add_thread( new boost::thread( boost::bind(permute,&ccdata,options.K,&keep,nperms_thread,gsl_rng_get(r),&permstats_thread[i]) ) );
	    }
	}
      tg.join_all();
      
      for(unsigned i=0;i<permstats_thread.size();++i)
	{
	  assert( permstats_thread[i].size() >= nperms_thread );
	  copy(permstats_thread[i].begin(),permstats_thread[i].end(),
	       std::back_inserter(permstats));
	}
      assert( permstats.size() == options.nperms );
      permstats_thread.clear();

    }
  sort(permstats.begin(),permstats.end());
  double perm_p = double( count_if(permstats.begin(),permstats.end(),
				   boost::bind(greater<double>(),_1,obs_esm ) ) )/double(options.nperms);
  double mean = gsl_stats_mean(&permstats[0],1,permstats.size());
  double sd = gsl_stats_sd(&permstats[0],1,permstats.size());
  double z = (obs_esm-mean)/sd;

  FILE * fp = fopen(options.outfile.c_str(),"a");
  int fd = fileno(fp);

  flock fd_lock = get_whole_flock();

  if (fcntl(fd,F_SETLKW,&fd_lock) == -1)
    {
      cerr << "ERROR: could not obtain lock on " << options.outfile << '\n';
      exit(10);
    }

  fprintf(fp,"%u\t%lf\t%lf\n",options.record_no,perm_p,z);

  fd_lock.l_type = F_UNLCK;
  if (fcntl(fd,F_UNLCK,&fd_lock) == -1)
    {
      cerr << "ERROR: could not release lock on " << options.outfile << '\n';
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

params parse_command_line(int argc, char ** argv)
{
  params rv;

  return rv;
}
