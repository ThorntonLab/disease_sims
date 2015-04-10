/*
  Helper fxns for cc2plink.R
 */

#include <Rcpp.h>
#include <vector>
#include <sstream>

using namespace Rcpp;
using namespace std;

//[[Rcpp::export(".reformatCCgenos")]]
DataFrame reformactCCgenos( const IntegerMatrix & m )
{
  vector< vector<char> > rv(2*m.ncol(),vector<char>(m.nrow(),'A'));

  for( int site = 0 ; site < m.ncol() ; ++site ) 
    {
      for( int ind = 0 ; ind < m.nrow() ; ++ind )
	{
	  int val = m(ind,site);
	  if(val == 0)
	    {
	      rv[2*site][ind]=rv[2*site+1][ind]='A';
	    }
	  else if (val == 1)
	    {
	      rv[2*site][ind]='A';
	      rv[2*site+1][ind]='T';
	    }
	  else if (val == 2)
	    {
	      rv[2*site][ind]=rv[2*site+1][ind]='T';
	    }
	  else 
	    {
	      stop("invalid value encountered");
	    }
	}
    }
  List temp(rv.size());
  Rcpp::CharacterVector colNames;
  for( int i = 0 ; i < temp.size() ; ++i) 
    {
      ostringstream NAME;
      NAME << 'm' << (i+1);
      colNames.push_back( NAME.str() );
      temp[i] = Rcpp::wrap(rv[i].begin(),rv[i].end());
    }
  temp.attr("names")=colNames;
  return Rcpp::DataFrame(temp);
}
