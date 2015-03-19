// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getRiskVariantMatrixDetails
Rcpp::List getRiskVariantMatrixDetails(const std::string& model, const std::string& popfile, const int64_t& popfile_offset, const unsigned& record_id_no, const double& dominance);
RcppExport SEXP diseaseSims_getRiskVariantMatrixDetails(SEXP modelSEXP, SEXP popfileSEXP, SEXP popfile_offsetSEXP, SEXP record_id_noSEXP, SEXP dominanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const std::string& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type popfile(popfileSEXP);
    Rcpp::traits::input_parameter< const int64_t& >::type popfile_offset(popfile_offsetSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type record_id_no(record_id_noSEXP);
    Rcpp::traits::input_parameter< const double& >::type dominance(dominanceSEXP);
    __result = Rcpp::wrap(getRiskVariantMatrixDetails(model, popfile, popfile_offset, record_id_no, dominance));
    return __result;
END_RCPP
}
// getRiskVariantMatrixDetails_Pheno
Rcpp::List getRiskVariantMatrixDetails_Pheno(const std::string& model, const std::string& popfile, const int64_t& popfile_offset, const std::string& phenofile, const int64_t& phenofile_offset, const unsigned& record_id_no, const double& dominance);
RcppExport SEXP diseaseSims_getRiskVariantMatrixDetails_Pheno(SEXP modelSEXP, SEXP popfileSEXP, SEXP popfile_offsetSEXP, SEXP phenofileSEXP, SEXP phenofile_offsetSEXP, SEXP record_id_noSEXP, SEXP dominanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const std::string& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type popfile(popfileSEXP);
    Rcpp::traits::input_parameter< const int64_t& >::type popfile_offset(popfile_offsetSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type phenofile(phenofileSEXP);
    Rcpp::traits::input_parameter< const int64_t& >::type phenofile_offset(phenofile_offsetSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type record_id_no(record_id_noSEXP);
    Rcpp::traits::input_parameter< const double& >::type dominance(dominanceSEXP);
    __result = Rcpp::wrap(getRiskVariantMatrixDetails_Pheno(model, popfile, popfile_offset, phenofile, phenofile_offset, record_id_no, dominance));
    return __result;
END_RCPP
}
// writeVpV1Data
void writeVpV1Data(const Rcpp::NumericMatrix& d, const std::string& outfilename, const unsigned& replicate_id, const bool& append);
RcppExport SEXP diseaseSims_writeVpV1Data(SEXP dSEXP, SEXP outfilenameSEXP, SEXP replicate_idSEXP, SEXP appendSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type outfilename(outfilenameSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type replicate_id(replicate_idSEXP);
    Rcpp::traits::input_parameter< const bool& >::type append(appendSEXP);
    writeVpV1Data(d, outfilename, replicate_id, append);
    return R_NilValue;
END_RCPP
}
// getEsizes
DataFrame getEsizes(const char * filename, const unsigned long& offset);
RcppExport SEXP diseaseSims_getEsizes(SEXP filenameSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const char * >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const unsigned long& >::type offset(offsetSEXP);
    __result = Rcpp::wrap(getEsizes(filename, offset));
    return __result;
END_RCPP
}
// sampleCCfromPop
List sampleCCfromPop(const char * popfilename, const unsigned long& offset, const char * phenofilename, const unsigned long& phenooffset, const unsigned& ncontrols, const unsigned& ncases, const double& case_proportion, const double& control_range, const unsigned& seed);
RcppExport SEXP diseaseSims_sampleCCfromPop(SEXP popfilenameSEXP, SEXP offsetSEXP, SEXP phenofilenameSEXP, SEXP phenooffsetSEXP, SEXP ncontrolsSEXP, SEXP ncasesSEXP, SEXP case_proportionSEXP, SEXP control_rangeSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const char * >::type popfilename(popfilenameSEXP);
    Rcpp::traits::input_parameter< const unsigned long& >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const char * >::type phenofilename(phenofilenameSEXP);
    Rcpp::traits::input_parameter< const unsigned long& >::type phenooffset(phenooffsetSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type ncontrols(ncontrolsSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type ncases(ncasesSEXP);
    Rcpp::traits::input_parameter< const double& >::type case_proportion(case_proportionSEXP);
    Rcpp::traits::input_parameter< const double& >::type control_range(control_rangeSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type seed(seedSEXP);
    __result = Rcpp::wrap(sampleCCfromPop(popfilename, offset, phenofilename, phenooffset, ncontrols, ncases, case_proportion, control_range, seed));
    return __result;
END_RCPP
}
// getCCblock
List getCCblock(const char * filename, const unsigned long& offset);
RcppExport SEXP diseaseSims_getCCblock(SEXP filenameSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const char * >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const unsigned long& >::type offset(offsetSEXP);
    __result = Rcpp::wrap(getCCblock(filename, offset));
    return __result;
END_RCPP
}
// getCCids
List getCCids(const char * filename, const unsigned long& offset);
RcppExport SEXP diseaseSims_getCCids(SEXP filenameSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const char * >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const unsigned long& >::type offset(offsetSEXP);
    __result = Rcpp::wrap(getCCids(filename, offset));
    return __result;
END_RCPP
}
// getPheno
NumericMatrix getPheno(const char * filename, const unsigned long& offset);
RcppExport SEXP diseaseSims_getPheno(SEXP filenameSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const char * >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const unsigned long& >::type offset(offsetSEXP);
    __result = Rcpp::wrap(getPheno(filename, offset));
    return __result;
END_RCPP
}
// writePVblock
void writePVblock(const char * outfilename, const char * indexfilename, const unsigned& recordno, DataFrame pvblock, const bool& append, const bool& gzip);
RcppExport SEXP diseaseSims_writePVblock(SEXP outfilenameSEXP, SEXP indexfilenameSEXP, SEXP recordnoSEXP, SEXP pvblockSEXP, SEXP appendSEXP, SEXP gzipSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const char * >::type outfilename(outfilenameSEXP);
    Rcpp::traits::input_parameter< const char * >::type indexfilename(indexfilenameSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type recordno(recordnoSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type pvblock(pvblockSEXP);
    Rcpp::traits::input_parameter< const bool& >::type append(appendSEXP);
    Rcpp::traits::input_parameter< const bool& >::type gzip(gzipSEXP);
    writePVblock(outfilename, indexfilename, recordno, pvblock, append, gzip);
    return R_NilValue;
END_RCPP
}
