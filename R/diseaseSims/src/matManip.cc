#include <matManip.hpp>
#include <algorithm>

using namespace std;
using namespace Rcpp;

/*
  Return a vector of integers:
  0 = column i is not a duplicate
  1 = column 1 is a duplicate of a least 1 column < 1

  Note: it would be tempting to return an Rcpp::LogicalVector
  here.  However, the Rcpp containers do not contain 
  reverse_iterator functions that are needed for the removal step 
  (see below).
 */
vector<int> columnsDuplicated( const IntegerMatrix & m )
{
  std::vector<int> rv(m.ncol(),0);
  unsigned nrow=m.nrow();
  IntegerMatrix::const_iterator i = m.begin();
  unsigned I = 0;
  for( ; i != m.end() ; i += nrow,++I )
    {
      unsigned J = I+1;
      for( IntegerMatrix::const_iterator j = i+nrow ;
	   j != m.end() ; j += nrow,++J )
	{
	  pair< IntegerMatrix::const_iterator,
		IntegerMatrix::const_iterator > pii = mismatch(i,i+nrow,
							       j);
	  if( pii.first == (i+nrow) && pii.second == j + nrow )
	    {
	      rv[J]=1;
	    }
	}
    }
  return rv;
}

/*
  Removes columns based on the 0/1 (false/true) data in lv.

  This works b/c an Rcpp matrix is "just a vector" with the dimensions
  added on top of it.

  We have to work backwards through rv so that we correctly remove the desired column.
*/
void removeDupColumns( Rcpp::IntegerMatrix & m,
		       const std::vector<int> & lv  )
{
  if(m.ncol() != lv.size())
    {
      Rcpp::stop("Error: number of columns in m must equal length of lv");
    }
   unsigned ncol = m.ncol(),nrow=m.nrow();

   unsigned removed = std::count(lv.begin(),lv.end(),1);

   for( std::vector<int>::const_reverse_iterator i = lv.rbegin() ;
	i != lv.rend() ; ++i )
     {
       if(*i==1)
	 {
	   std::vector<int>::iterator::difference_type d = std::distance(lv.begin(),i.base());
	   m.erase( m.begin() + (d-1)*nrow, m.begin() + d*nrow );
	 }
     }

   m.attr("dim") = Dimension(nrow,ncol-removed);
}
