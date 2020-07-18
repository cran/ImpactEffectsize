// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// c_gmd
double c_gmd(Rcpp::NumericVector j);
RcppExport SEXP _ImpactEffectsize_c_gmd(SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(c_gmd(j));
    return rcpp_result_gen;
END_RCPP
}
// c_median
double c_median(Rcpp::NumericVector x);
RcppExport SEXP _ImpactEffectsize_c_median(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_median(x));
    return rcpp_result_gen;
END_RCPP
}
// c_quantile
NumericVector c_quantile(NumericVector x, NumericVector probs);
RcppExport SEXP _ImpactEffectsize_c_quantile(SEXP xSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_quantile(x, probs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ImpactEffectsize_c_gmd", (DL_FUNC) &_ImpactEffectsize_c_gmd, 1},
    {"_ImpactEffectsize_c_median", (DL_FUNC) &_ImpactEffectsize_c_median, 1},
    {"_ImpactEffectsize_c_quantile", (DL_FUNC) &_ImpactEffectsize_c_quantile, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ImpactEffectsize(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}