// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// C_get_from_pos
List C_get_from_pos(Environment m, IntegerMatrix pos_k);
RcppExport SEXP dyngame_C_get_from_pos(SEXP mSEXP, SEXP pos_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type m(mSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type pos_k(pos_kSEXP);
    rcpp_result_gen = Rcpp::wrap(C_get_from_pos(m, pos_k));
    return rcpp_result_gen;
END_RCPP
}
// C_RowSums_VectorList
NumericVector C_RowSums_VectorList(NumericVector vec, IntegerVector ncols);
RcppExport SEXP dyngame_C_RowSums_VectorList(SEXP vecSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(C_RowSums_VectorList(vec, ncols));
    return rcpp_result_gen;
END_RCPP
}
// C_RowMaxs_VectorList
NumericVector C_RowMaxs_VectorList(NumericVector vec, IntegerVector ncols);
RcppExport SEXP dyngame_C_RowMaxs_VectorList(SEXP vecSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(C_RowMaxs_VectorList(vec, ncols));
    return rcpp_result_gen;
END_RCPP
}
// C_which_RowMaxs_VectorList
IntegerVector C_which_RowMaxs_VectorList(NumericVector vec, IntegerVector ncols);
RcppExport SEXP dyngame_C_which_RowMaxs_VectorList(SEXP vecSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(C_which_RowMaxs_VectorList(vec, ncols));
    return rcpp_result_gen;
END_RCPP
}
