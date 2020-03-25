#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

class myclass {
private:
  arma::vec x;
  std::vector< arma::mat * > y;
public:
  myclass(const NumericVector & dat, const ListOf< NumericMatrix > & LM) : 
  x((double *) &dat[0], dat.size(), false, true)
  {
    for (int i = 0; i < LM.size(); ++i) {
      y.push_back(new arma::mat((double *) &LM[i][0],  LM[i].nrow(), LM[i].ncol(), false, true));
    }
    
    return;
  }
  
  std::vector< double > mean() const {
    std::vector< double > ans(1 + y.size());
    ans[0] = arma::mean(x);
    for (int i = 0; i < y.size(); ++i)
      ans[i + 1] = arma::mean(arma::mean((*y[i])));
      
    return ans;
  }
  
  std::vector< double > sd() const {
    std::vector< double > ans(1 + y.size());
    ans[0] = arma::stddev(x);
    for (int i = 0; i < y.size(); ++i)
      ans[i + 1] = arma::stddev(arma::mean((*y[i])));
    
    return ans;
  }
};

// [[Rcpp::export]]
SEXP mywrapper(const NumericVector & x, const ListOf< NumericMatrix > & LM) {
  
  Rcpp::XPtr< myclass > x_ptr(
    new myclass(x, LM)
  );
  
  return x_ptr;
   
}


// [[Rcpp::export]]
NumericVector wrapmean(SEXP x) {
  Rcpp::XPtr< myclass > p(x);
  return Rcpp::wrap(p->mean());
}

// [[Rcpp::export]]
NumericVector wrapsd(SEXP x) {
  Rcpp::XPtr< myclass > p(x);
  return Rcpp::wrap(p->sd());
}

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
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
dat <- runif(1e5)
M   <- replicate(10, matrix(runif(1e8/10), ncol = 20), simplify = FALSE)
ptr <- mywrapper(dat, M)
ans0_mean <- wrapmean(ptr)
ans1_mean <- c(mean(dat), sapply(M, mean))
ans0_mean - ans1_mean

ans0_sd <- wrapsd(ptr)
ans1_sd <- c(sd(dat), sapply(M, function(i) sd(colMeans(i))))
ans0_sd - ans1_sd
*/
