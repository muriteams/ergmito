#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double exp2(double x) {
  x = 1.0 + x / 1024.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x;
  return x;
}

// [[Rcpp::export]]
double exp1(double x) {
  return exp(x);
}

/***R
plot(microbenchmark(
  exp = exp1(5),
  exp1024 = exp2(5), times = 10000
))
*/