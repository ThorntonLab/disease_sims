#include <matManip.hpp>
#include <algorithm>

using namespace std;
using namespace Rcpp;

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
