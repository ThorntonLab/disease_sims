#include <Rcpp.h>
#include <vector>
#include <utility>
#include <string>
#include <cstdint>
#include <zlib.h>
#include <fwdpp/IO.hpp>
#include <diseaseSims/util.hpp>
#include <diseaseSims/traitValues.hpp>

using namespace std;

// Details of how to get a genotype matrix for risk variants
//[[Rcpp::export(".getRiskVariantMatrixDetails")]]
Rcpp::List getRiskVariantMatrixDetails( const std::string & model,
					const std::string & popfile,
					const int64_t & popfile_offset,
					const unsigned & record_id_no,
					const double & dominance = 0.)
{
  gzFile gzin = gzopen(popfile.c_str(),"rb");
  if( gzin == NULL ) 
    {
      Rcpp::Rcerr << "Error, " << popfile
		  << " could not be opened for reading.\n";
      return Rcpp::List();
    }

  gzseek(gzin,popfile_offset, SEEK_SET);

  Gfxn_t dipG = std::bind(TFL2013g(),std::placeholders::_1,std::placeholders::_2);
  if( model == "additive" )
    {
      dipG = std::bind(additiveg(),std::placeholders::_1,std::placeholders::_2);
    }
  else if (model == "multipicative")
    {
      dipG = std::bind(multiplicative_phenotype(),std::placeholders::_1,std::placeholders::_2);
    }
  else if (model == "popgen")
    {
      dipG = std::bind(popgen_phenotype(),std::placeholders::_1,std::placeholders::_2,dominance);
    }

  popstruct pop = readPop(gzin);
  // mlist mutations;
  // glist gametes;
  // dipvector diploids;
  // KTfwd::read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
  gzclose(gzin);

  unsigned RISKMUTIDX=0;
  vector<pair<mlist::iterator,unsigned> > risk_indexes;
  for( auto i = pop.mutations.begin();i!=pop.mutations.end();++i )
    {
      if( ! i->neutral )
	{
	  risk_indexes.push_back( make_pair(i,RISKMUTIDX++) );
	}
    }
  auto Gvals = getG(pop.diploids,dipG);

  Rcpp::IntegerMatrix genos(pop.diploids.size(),RISKMUTIDX);

  for( unsigned ind = 0 ; ind < pop.diploids.size() ; ++ind )
    {
      vmcount_t vmc = get_mut_counts(pop.diploids[ind].first,pop.diploids[ind].second);
      for( unsigned i = 0 ; i < vmc.size() ; ++i )
	{
	  auto __itr = find_if( risk_indexes.begin(), risk_indexes.end(),[&vmc,&i](const pair<mlist::iterator,unsigned> & __p) {
	      return __p.first == vmc[i].first;
	    });
	  genos(ind,__itr->second) += vmc[i].second;
	}
    }
  return Rcpp::List::create(Rcpp::Named("G") = Gvals,
			    Rcpp::Named("genos") = genos);
}

