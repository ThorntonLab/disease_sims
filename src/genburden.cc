/*
 Read the population file from TFL2013_ind.cc and calculate some burden summary stats
 User must supply:
 #1) The number of replicates to process
 #2) The name of the population file
 #3) The name of the desired output file
*/


#include <iostream>
#include <fwdpp/diploid.hh>
#include <utility>
#include <iostream>
#include <fstream>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/program_options.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <simindex.hpp>

#include <zlib.h>

#include <diseaseSims/mutation_with_age.hpp>
#include <fwdpp/sugar/serialization.hpp>
using namespace std;
using namespace boost::program_options;
using namespace boost::interprocess;
using namespace boost::accumulators;
using namespace KTfwd;

struct popparams 
{
  unsigned reps;
  string popfile,ofile; 
  popparams(void);
};

popparams::popparams(void): reps(250),
			    popfile(string()),
			    ofile(string())
{
}

popparams parse_command_line(const int & argc,
			    char **argv);

popparams parse_command_line(const int & argc, char ** argv)
{
  popparams rv;

  options_description desc("Summarize genetic burden from TFL2013");

  desc.add_options()
    ("help,h", "Produce help message")
    ("reps,r",value<unsigned>(&rv.reps)->default_value(250),"Number of replicates.")
    ("popfile,p",value<string>(&rv.popfile)->default_value(string()),"Population file output by TFL2013_ind")
    ("ofile,o",value<string>(&rv.ofile)->default_value(string()),"Name of output file")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);


  if(argc == 1 || vm.count("help"))
    {
      cerr << desc << '\n';
      exit(0);
    }

  if(! vm.count("popfile") )
    {
      cerr << "Error, no popfile specified\n";
      exit(EXIT_FAILURE);
    }

  if(! vm.count("ofile") )
    {
      cerr << "Error, no output file specified\n";
      exit(EXIT_FAILURE);
    }
  return rv;
};


int main(int argc, char **argv )
{
  
  popparams params = parse_command_line(argc,argv);  
  gzFile gzin = gzopen(params.popfile.c_str(),"rb");
  ofstream output;
  output.open(params.ofile.c_str());
  output << "mat_mean_nmut"<<' '<< "mat_var_nmut"<<' '<<"mat_skew_nmut"<<' '<<"mat_kurt_nmut" <<' '
	 <<"pat_mean_nmut"<<' '<<"pat_var_nmut"<<' '<<"pat_skew_nmut"<<' '<<"pat_kurt_nmut" <<' '
	 << "mat_mean_eff"<<' '<< "mat_var_eff"<<' '<<"mat_skew_eff"<<' '<<"mat_kurt_eff" <<' '
	 <<"pat_mean_eff"<<' '<<"pat_var_eff" <<' '<<"pat_skew_eff"<<' '<<"pat_kurt_eff" <<'\n';
  for ( unsigned i = 0; i <params.reps; ++i ) {
    poptype pop(0);
    KTfwd::gzdeserialize()(gzin,pop,KTfwd::mutation_reader<TFLmtype>());
    //read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
    vector<unsigned> nmom(pop.diploids.size());
    vector<unsigned> ndad(pop.diploids.size());
    vector<double> emom(pop.diploids.size());
    vector<double> edad(pop.diploids.size());
    std::size_t j=0;
    for(const auto & dip : pop.diploids)
      {
	nmom[j]=pop.gametes[dip.first].smutations.size();
	ndad[j]=pop.gametes[dip.second].smutations.size();
	emom[j]={0.},edad[j]={0.};
	for(const auto m : pop.gametes[dip.first].smutations) emom[j]+= pop.mutations[m].s;
	for(const auto m : pop.gametes[dip.second].smutations) edad[j]+= pop.mutations[m].s;
	++j;
      }
    accumulator_set<double, stats<tag::mean, tag::variance, tag::skewness, tag::kurtosis > > nmacc;
    for_each(nmom.begin(),nmom.end(),bind<void>(ref(nmacc),_1));
    
    accumulator_set<double, stats<tag::mean, tag::variance, tag::skewness, tag::kurtosis > > ndacc;
    for_each(ndad.begin(),ndad.end(),bind<void>(ref(ndacc),_1));

    accumulator_set<double, stats<tag::mean, tag::variance, tag::skewness, tag::kurtosis > > emacc;
    for_each(emom.begin(),emom.end(),bind<void>(ref(emacc),_1));

    accumulator_set<double, stats<tag::mean, tag::variance, tag::skewness, tag::kurtosis > > edacc;
    for_each(edad.begin(),edad.end(),bind<void>(ref(edacc),_1));

    output << mean( nmacc ) <<' '<< variance(nmacc) <<' '<< skewness(nmacc) <<' '<< kurtosis(nmacc) <<' '
	   << mean( nmacc ) <<' '<< variance(ndacc) <<' '<< skewness(ndacc) <<' '<< kurtosis(ndacc) <<' '
	   << mean( emacc ) <<' '<< variance(emacc) <<' '<< skewness(emacc) <<' '<< kurtosis(emacc) <<' '
	   << mean( edacc ) <<' '<< variance(edacc) <<' '<< skewness(edacc) <<' '<< kurtosis(edacc) <<'\n';
    
  }
  gzclose(gzin);
  output.close();
}

