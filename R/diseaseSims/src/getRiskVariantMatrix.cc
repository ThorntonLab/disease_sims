#include <Rcpp.h>
#include <vector>
#include <utility>
#include <string>
#include <cstdint>
#include <algorithm>
#include <sstream>
#include <set>
#include <zlib.h>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <fwdpp/IO.hpp>
#include <diseaseSims/util.hpp>
#include <diseaseSims/traitValues.hpp>

using namespace std;

Gfxn_t setModel( const std::string & model, const double & dominance )
{
  Gfxn_t dipG = std::bind(TFL2013g(),std::placeholders::_1,std::placeholders::_2);
  if( model == "additive" )
    {
      dipG = std::bind(additiveg(),std::placeholders::_1,std::placeholders::_2);
    }
  else if (model == "multiplicative")
    {
      dipG = std::bind(multiplicative_phenotype(),std::placeholders::_1,std::placeholders::_2);
    }
  else if (model == "popgen")
    {
      dipG = std::bind(popgen_phenotype(),std::placeholders::_1,std::placeholders::_2,dominance);
    }
  return dipG;
}

/*
  Return value is sorted by:
  1. Decreasing mutation frequency
  2. Within frequency class, by decreasing abs(effect size)
 */
vector<pair<mlist::iterator,unsigned> > getRiskIndexes( mlist & mutations )
{
  vector<mlist::iterator> v;
  set<unsigned,greater<unsigned> > counts;
  for( auto i = mutations.begin();i!=mutations.end();++i ) 
    {
      if(!i->neutral)
	{
	  v.push_back(i);
	}
    }
  std::sort(v.begin(),v.end(),[&counts]( mlist::iterator const & i, mlist::iterator const & j ) { 
      counts.insert(i->n);
      counts.insert(j->n);
      return i->n > j->n; 
    } );
  //Sort by decreasing |effect size| w/in each freq class
  for_each( counts.begin(), counts.end(),
	    [&v]( const unsigned & __u ) {
	      //First element at this freq.
	      auto __beg = find_if(v.begin(),v.end(),[&__u](const mlist::iterator & __i) { return __i->n == __u; });
	      //Last element at this freq.
	      auto __end = find_if(v.rbegin(),v.rend(),[&__u](const mlist::iterator & __i) { return __i->n == __u; });
	      sort(__beg,__end.base(),[](mlist::iterator const & __i, mlist::iterator const & __j) {
		  return abs(__i->s) > abs(__j->s);
		});
	    } );
  vector<pair<mlist::iterator,unsigned> > risk_indexes;
  unsigned RISKMUTIDX=0;
  for( auto i = v.begin() ; i != v.end() ; ++i )
    {
      risk_indexes.push_back( make_pair(*i,RISKMUTIDX++) );
    }
  return risk_indexes;
}

vector<double> getEsizes( const vector<pair<mlist::iterator,unsigned> > & risk_indexes )
{
  vector<double> rv;
  for_each(risk_indexes.begin(),risk_indexes.end(),[&rv](const pair<mlist::iterator,unsigned> & __p) {
      rv.push_back(__p.first->s);
    });
  return rv;
}

std::vector<std::int8_t> columnsDuplicated( const std::vector<std::vector<unsigned> > & v )
{
  std::vector<std::int8_t> rv(v.size(),0);
  unsigned I = 0;
  for( ; I < rv.size()-1 ; ++I )
    {
      for( unsigned J = I + 1 ; J < rv.size() ; ++J )
	{
	  auto XX = std::mismatch( v[I].begin(),v[I].end(), v[J].begin() );
	  if( XX.first == v[I].end() && XX.second == v[J].end() )
	    {
	      rv[J]=1;
	    }
	}
    }
  return rv;
}

Rcpp::DataFrame MakeRiskMatrix( const dipvector & diploids,
				const std::vector<double> & trait_vals,
				const vector<pair<mlist::iterator,unsigned> > & risk_indexes )
/*
  Problem: R/Rcpp matrices may not be able to hold enough data: http://stackoverflow.com/questions/9984283/maximum-size-of-a-matrix-in-r
  Solution: Romain's advice from:
  http://stackoverflow.com/questions/23865210/how-to-convert-stdvectorstdvectordouble-to-rcppdataframe-or-rcppnume

  Romain implies that it may be important to set the dimnames.  It is, otherwise things go south with the return value.
*/
{
  std::vector<std::vector<unsigned> > temp(risk_indexes.size(),
					   std::vector<unsigned>(diploids.size(),0u));
  for( unsigned ind = 0 ; ind < diploids.size() ; ++ind )
    {
      vmcount_t vmc = get_mut_counts(diploids[ind].first,diploids[ind].second);
      for( unsigned i = 0 ; i < vmc.size() ; ++i )
	{
	  auto __itr = find_if( risk_indexes.begin(), risk_indexes.end(),[&vmc,&i](const pair<mlist::iterator,unsigned> & __p) {
	      return __p.first == vmc[i].first;
	    });
	  temp[__itr->second][ind] = vmc[i].second;
	}
    }
  /* 
     Let's remove identical columns.
     lm() and the like will do that anyways.
     Also, pruning here lowers total data size, etc.
  */
  /*
  auto dup = columnsDuplicated(temp);
  auto nremoved = count(dup.begin(),dup.end(),1);
  if( nremoved ) //any dups??
    {
      for( auto itr = dup.rbegin() ; itr != dup.rend() ; ++itr )
	{
	  if(*itr)
	    {
	      auto d = std::distance(dup.begin(),itr.base());
	      temp.erase(temp.begin() + d-1);
	      esizes.erase(esizes.begin() + d-1);
	    }
	}
    }
  */
  Rcpp::List temp2(temp.size()+1);
  temp2[0] = Rcpp::wrap( trait_vals.begin(),trait_vals.end() );
  Rcpp::CharacterVector colNames;
  colNames.push_back("trait");
  unsigned i = 0;
  for( ; i < temp.size() ; ++i) 
    {
      ostringstream NAME;
      NAME << 'V' << (i+1);
      colNames.push_back( NAME.str() );
      temp2[i+1] = Rcpp::wrap(temp[i].begin(),temp[i].end());
    }
  temp2.attr("names")=colNames;
  return Rcpp::DataFrame(temp2);
}

// Details of how to get a genotype matrix for risk variants
//[[Rcpp::export(".getRiskVariantMatrixDetails")]]
Rcpp::List getRiskVariantMatrixDetails( const std::string & model,
					const std::string & popfile,
					const int64_t & popfile_offset,
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

  Gfxn_t dipG = setModel(model,dominance);

  popstruct pop = readPop(gzin);
  gzclose(gzin);

  vector<pair<mlist::iterator,unsigned> > risk_indexes = getRiskIndexes(pop.mutations);
  auto Gvals = getG(pop.diploids,dipG);
  auto genos = MakeRiskMatrix(pop.diploids,Gvals,risk_indexes);
  return Rcpp::List::create(Rcpp::Named("esizes") = getEsizes(risk_indexes),
			    Rcpp::Named("genos") = genos);
 }
				    
// Details of how to get a genotype matrix for risk variants
//[[Rcpp::export(".getRiskVariantMatrixDetailsPheno")]]
Rcpp::List getRiskVariantMatrixDetails_Pheno( const std::string & model,
					      const std::string & popfile,
					      const int64_t & popfile_offset,
					      const std::string & phenofile,
					      const int64_t & phenofile_offset,
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

  Gfxn_t dipG = setModel(model,dominance);

  popstruct pop = readPop(gzin);
  gzclose(gzin);

  gzin = gzopen(phenofile.c_str(),"rb");
  gzseek(gzin,phenofile_offset,SEEK_SET);
  unsigned ndips = 0;
  gzread(gzin,&ndips,sizeof(unsigned));
  vector<double> phenotypes(ndips);
  for( unsigned i = 0 ; i < ndips ; ++i )
    {
      double G,E;
      gzread(gzin,&G,sizeof(double));
      gzread(gzin,&E,sizeof(double));
      phenotypes[i] = G + E;
    }
  gzclose(gzin);
  vector<pair<mlist::iterator,unsigned> > risk_indexes = getRiskIndexes(pop.mutations);
  auto genos = MakeRiskMatrix(pop.diploids,phenotypes,risk_indexes);
  return Rcpp::List::create(Rcpp::Named("esizes") = getEsizes(risk_indexes),
			    Rcpp::Named("genos") = genos);
 }
