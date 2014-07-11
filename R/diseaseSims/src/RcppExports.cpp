// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getEsizes
DataFrame getEsizes(const char * filename, const unsigned long& offset);
RcppExport SEXP diseaseSims_getEsizes(SEXP filenameSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const char * >::type filename(filenameSEXP );
        Rcpp::traits::input_parameter< const unsigned long& >::type offset(offsetSEXP );
        DataFrame __result = getEsizes(filename, offset);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// getCCblock
List getCCblock(const char * filename, const unsigned long& offset);
RcppExport SEXP diseaseSims_getCCblock(SEXP filenameSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const char * >::type filename(filenameSEXP );
        Rcpp::traits::input_parameter< const unsigned long& >::type offset(offsetSEXP );
        List __result = getCCblock(filename, offset);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// getPheno
NumericMatrix getPheno(const char * filename, const unsigned long& offset);
RcppExport SEXP diseaseSims_getPheno(SEXP filenameSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const char * >::type filename(filenameSEXP );
        Rcpp::traits::input_parameter< const unsigned long& >::type offset(offsetSEXP );
        NumericMatrix __result = getPheno(filename, offset);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// writePVblock
void writePVblock(const char * outfilename, const char * indexfilename, const unsigned& recordno, DataFrame pvblock, const bool& append = true, const bool& gzip = false);
RcppExport SEXP diseaseSims_writePVblock(SEXP outfilenameSEXP, SEXP indexfilenameSEXP, SEXP recordnoSEXP, SEXP pvblockSEXP, SEXP appendSEXP, SEXP gzipSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const char * >::type outfilename(outfilenameSEXP );
        Rcpp::traits::input_parameter< const char * >::type indexfilename(indexfilenameSEXP );
        Rcpp::traits::input_parameter< const unsigned& >::type recordno(recordnoSEXP );
        Rcpp::traits::input_parameter< DataFrame >::type pvblock(pvblockSEXP );
        Rcpp::traits::input_parameter< const bool& >::type append(appendSEXP );
        Rcpp::traits::input_parameter< const bool& >::type gzip(gzipSEXP );
        writePVblock(outfilename, indexfilename, recordno, pvblock, append, gzip);
    }
    return R_NilValue;
END_RCPP
}