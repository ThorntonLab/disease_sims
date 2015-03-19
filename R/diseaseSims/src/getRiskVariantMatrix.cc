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

Gfxn_t setModel( const std::string & model, const double & dominance )
{
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
  return dipG;
}

vector<pair<mlist::iterator,unsigned> > getRiskIndexes( mlist & mutations )
{
  vector<pair<mlist::iterator,unsigned> > risk_indexes;
  unsigned RISKMUTIDX=0;
  for( auto i = mutations.begin();i!=mutations.end();++i )
    {
      if( ! i->neutral )
	{
	  risk_indexes.push_back( make_pair(i,RISKMUTIDX++) );
	}
    }
  return risk_indexes;
}

Rcpp::IntegerMatrix MakeRiskMatrik( const dipvector & diploids,
				    const vector<pair<mlist::iterator,unsigned> > & risk_indexes )
{
  Rcpp::IntegerMatrix genos(diploids.size(),risk_indexes.size());
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
  return genos;
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

  Rcpp::IntegerMatrix genos = MakeRiskMatrik(pop.diploids,risk_indexes);

  return Rcpp::List::create(Rcpp::Named("trait") = Gvals,
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

  Rcpp::IntegerMatrix genos = MakeRiskMatrik(pop.diploids,risk_indexes);

  return Rcpp::List::create(Rcpp::Named("trait") = phenotypes,
			    Rcpp::Named("genos") = genos);
 }

//[[Rcpp::export(".writeVpV1Data")]]
 void writeVpV1Data( const Rcpp::NumericMatrix & d,
		     const std::string & outfilename,
		     const unsigned & replicate_id,
		     const bool & append ) 
 {
   std::string __append = (append) ? "a" : "w";
   gzFile gzout = gzopen(outfilename.c_str(),__append.c_str());

   for(int row = 0 ; row < d.nrow() ; ++row )
     {
       gzprintf(gzout,"%u\t%lf\t%lf\t%lf\n",replicate_id,d(row,0),d(row,1),d(row,2));
     }
   gzclose(gzout);
 }
