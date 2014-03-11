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
using namespace boost::program_options;

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
  /*
    string ccindexfile,ccfile,outfile;
    unsigned record_no,K,max_controls,max_cases,nperms,nthreads,seed;
    double minfreq,maxfreq,rsq_cutoff;
  */
  options_description desc("Calculate the \"ESM\" statistic of Thornton, Foran and Long(2013) PLoS Genetics 9(2): e1003258.\nThe pvalue based on nperms permutations of case/control status is output,\nalong with the z-score corresponding to the statistic,\nfollowing the approach described in the PLoS Genetics paper (read it for details).");
  desc.add_options()
    ("help,h", "Produce help message")
    ("ccindex,c",value<string>(&rv.ccindexfile),"Name of index file generated by make_case_control.")
    ("ccfile,C",value<string>(&rv.ccfile),"Name of output file generated by make_case_control.")
    ("outfile,o",value<string>(&rv.outfile),"Name of output file.")
    ("recordno,r",value<unsigned>(&rv.record_no),"Record number to look up in ccindex.")
    ("maxmarkers,K",value<unsigned>(&rv.K),"Max number of SNPs to use for ESM statistic.")
    ("maxcontrols",value<unsigned>(&rv.max_controls)->default_value(numeric_limits<unsigned>::max()),"Max number of control individuals to use in analysis")
    ("maxcases",value<unsigned>(&rv.max_cases)->default_value(numeric_limits<unsigned>::max()),"Max number of case individuals to use in analysis")
    ("nperms,n",value<unsigned>(&rv.nperms)->default_value(10000),"Number of permutations to perform")
    ("nthreads,t",value<unsigned>(&rv.nthreads)->default_value(1),"Number of permutations to perform")
    ("seed,s",value<unsigned>(&rv.seed)->default_value(0),"Random number seed")
    ("minfreq,m",value<double>(&rv.minfreq)->default_value(0.),"Throw out a marker unless its minor allele frequency in controls is >= minfreq")
    ("maxfreq,m",value<double>(&rv.maxfreq)->default_value(0.05),"Throw out a marker unless its minor allele frequency in controls is < maxfreq")
    ("rsq-cutoff,R",value<double>(&rv.rsq_cutoff)->default_value(0.8),"The program calculates r^2 for all pairs of sites (based on diploid genotypes).  If r^2 >= rsq-cutoff, then only one of the two sites in a pair is retained for calculation of the ESM statistic")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if(argc==1 || vm.count("help"))
    {
      cout << desc << '\n';
      exit(0);
    }

  bool options_fail=false;
  if (! vm.count("ccindex") )
    {
      options_fail=true;
      cerr << "Error: --ccindex/-c option must be specified with an argument to specify the case/control index file\n";
    }
  if (! vm.count("ccfile") )
    {
      options_fail=true;
      cerr << "Error: --ccfile/-C option must be specified with an argument to specify the case/control data file\n";
    }
  if (! vm.count("outfile"))
    {
      options_fail = true;
      cerr << "Error: --outfile/-o must be used to specify outfile\n";
    }
  if (! vm.count("recordno"))
    {
      options_fail = true;
      cerr << "Error: --recordno/-r must be used to specify a record ID to looup in the case/control index file\n";
    }

  if(options_fail)
    {
      cerr << "Please use --help/-h to see complete list of program options\n";
      exit(10);
    }

  cerr <<"here\n";
  exit(0);			   
  return rv;
}
