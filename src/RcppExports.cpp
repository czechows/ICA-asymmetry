// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>
// [[Rcpp::depends(dlib)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// ICA
RcppExport SEXP ICA(const NumericMatrix& XX, NumericVector& mm, NumericMatrix& WW, double& c, int gauss_noise, bool generalized, double minimum, double accuracy);
RcppExport SEXP ICAA_ICA(SEXP XXSEXP, SEXP mmSEXP, SEXP WWSEXP, SEXP cSEXP, SEXP gauss_noiseSEXP, SEXP generalizedSEXP, SEXP minimumSEXP, SEXP accuracySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type WW(WWSEXP);
    Rcpp::traits::input_parameter< double& >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type gauss_noise(gauss_noiseSEXP);
    Rcpp::traits::input_parameter< bool >::type generalized(generalizedSEXP);
    Rcpp::traits::input_parameter< double >::type minimum(minimumSEXP);
    Rcpp::traits::input_parameter< double >::type accuracy(accuracySEXP);
    rcpp_result_gen = Rcpp::wrap(ICA(XX, mm, WW, c, gauss_noise, generalized, minimum, accuracy));
    return rcpp_result_gen;
END_RCPP
}
