// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>
// [[Rcpp::depends(dlib)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// ICA
RcppExport SEXP ICA(const NumericMatrix& XX, NumericVector& mm, NumericMatrix& WW);
RcppExport SEXP ICAA_ICA(SEXP XXSEXP, SEXP mmSEXP, SEXP WWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type WW(WWSEXP);
    rcpp_result_gen = Rcpp::wrap(ICA(XX, mm, WW));
    return rcpp_result_gen;
END_RCPP
}
