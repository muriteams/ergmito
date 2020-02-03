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
std::vector< std::vector< int > > coordinates(const std::vector< int > & x, int n) {
  
  std::vector< std::vector< int > > res(n);
  int m = x.size(), xmax = n*n - 1;
  
  for (int i = 0; i < m; ++i) {
    
    if (x[i] > xmax | x[i] < 0)
      stop("Outside range.");
    
    res[x[i] % n].push_back(x[i] / n);
    
  }
  
  return res;
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
coordinates(c(0, 3, 6), 3)
coordinates(c(1, 4, 7), 3)
coordinates(c(2, 5, 8), 3)
*/
