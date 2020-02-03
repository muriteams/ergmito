#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
bool list2map(List x) {
  
  return x.containsElementNamed("Casa");

  
}

/***R

list2map(list(Casa=1))

*/