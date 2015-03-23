#ifndef __MATMANIP_HPP__
#define __MATMANIP_HPP__

#include <vector>
#include <Rcpp.h>

std::vector<int> columnsDuplicated( const Rcpp::IntegerMatrix & m );

void removeColumns( Rcpp::IntegerMatrix & m,
		    const std::vector<int> & lv  );

#endif
