#include <Rcpp.h>
#include <vector>
#include <utility>
#include <string>
#include <cstdint>
#include <zlib.h>
#include <fwdpp/IO.hpp>
#include <diseaseSims/traitValues.hpp>

using namespace std;

using Gfxn_t = std::function<double(const glist::const_iterator &,
				    const glist::const_iterator &)>;
using diploid_t = std::pair<glist::const_iterator,glist::const_iterator>;
using vmcount_t = vector<pair<mlist::iterator,int8_t> >;
vector<double> getG( const dipvector & diploids,
		     const Gfxn_t & dipG )
{
  vector<double> rv;
  for_each( diploids.begin(),diploids.end(),[&rv,&dipG](const diploid_t & __d ) { rv.push_back( dipG(__d.first,__d.second) ); } );
  return rv;
}

vmcount_t get_mut_counts( const glist::const_iterator & g1,
			  const glist::const_iterator & g2 );

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

  mlist mutations;
  glist gametes;
  dipvector diploids;
  KTfwd::read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
  gzclose(gzin);

  unsigned RISKMUTIDX=0;
  vector<pair<mlist::iterator,unsigned> > risk_indexes;
  for( auto i = mutations.begin();i!=mutations.end();++i )
    {
      if( ! i->neutral )
	{
	  risk_indexes.push_back( make_pair(i,RISKMUTIDX++) );
	}
    }
  auto Gvals = getG(diploids,dipG);

  Rcpp::IntegerMatrix genos(diploids.size(),RISKMUTIDX);

  for( unsigned ind = 0 ; ind < diploids.size() ; ++ind )
    {
      vmcount_t vmc = get_mut_counts(diploids[ind].first,diploids[ind].second);
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

vmcount_t get_mut_counts( const glist::const_iterator & g1,
			  const glist::const_iterator & g2 )
{
  vmcount_t rv;

  auto updater = [&rv](const mlist::iterator & __mut) {
    auto __itr =  find_if(rv.begin(),rv.end(),[&__mut](const pair<mlist::iterator,unsigned> & __p) {
	return __p.first == __mut;
      } );
    if(__itr == rv.end())
      rv.push_back(make_pair(__mut,1u));
    else
      __itr->second++;
  };
  for_each( g1->smutations.begin(), g1->smutations.end(),updater );
  for_each( g2->smutations.begin(), g2->smutations.end(),updater );
  return rv;
}
