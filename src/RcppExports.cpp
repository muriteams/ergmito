// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "lergm_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// exact_loglik
arma::vec exact_loglik(const arma::mat& x, const arma::rowvec& params, const arma::rowvec& weights, const arma::mat& statmat);
RcppExport SEXP _lergm_exact_loglik(SEXP xSEXP, SEXP paramsSEXP, SEXP weightsSEXP, SEXP statmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type statmat(statmatSEXP);
    rcpp_result_gen = Rcpp::wrap(exact_loglik(x, params, weights, statmat));
    return rcpp_result_gen;
END_RCPP
}
// make_sets
vecint make_sets(int n);
RcppExport SEXP _lergm_make_sets(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(make_sets(n));
    return rcpp_result_gen;
END_RCPP
}
// powerset
SEXP powerset(int n, bool force);
RcppExport SEXP _lergm_powerset(SEXP nSEXP, SEXP forceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type force(forceSEXP);
    rcpp_result_gen = Rcpp::wrap(powerset(n, force));
    return rcpp_result_gen;
END_RCPP
}
// print_powerset
int print_powerset(SEXP sets);
RcppExport SEXP _lergm_print_powerset(SEXP setsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sets(setsSEXP);
    rcpp_result_gen = Rcpp::wrap(print_powerset(sets));
    return rcpp_result_gen;
END_RCPP
}
// wrap_powerset
List wrap_powerset(SEXP sets, int from, int to, int n);
RcppExport SEXP _lergm_wrap_powerset(SEXP setsSEXP, SEXP fromSEXP, SEXP toSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sets(setsSEXP);
    Rcpp::traits::input_parameter< int >::type from(fromSEXP);
    Rcpp::traits::input_parameter< int >::type to(toSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(wrap_powerset(sets, from, to, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lergm_exact_loglik", (DL_FUNC) &_lergm_exact_loglik, 4},
    {"_lergm_make_sets", (DL_FUNC) &_lergm_make_sets, 1},
    {"_lergm_powerset", (DL_FUNC) &_lergm_powerset, 2},
    {"_lergm_print_powerset", (DL_FUNC) &_lergm_print_powerset, 1},
    {"_lergm_wrap_powerset", (DL_FUNC) &_lergm_wrap_powerset, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_lergm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
