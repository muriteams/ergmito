#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(name = "geodesic_cpp")]]
IntegerMatrix geodesic(const arma::imat & x) {
  
  arma::imat res(x.n_rows, x.n_cols);
  res.fill(-1);
  arma::imat res_tmp = x;
  
  // List of indices to check
  int n = x.n_rows;
  unsigned int changes = 0u;
  for (int iter = 0; iter < n; ++iter) {
    changes = 0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        if (i != j && res_tmp.at(i, j) != 0 && (res.at(i, j) == -1)) 
          res.at(i, j) = iter + 1, ++changes;
      }
      
    // // No more changes to record
    // if (changes == 0u)
    //   break;
    
    // Powermatrix
    res_tmp *= x;
  }
  // print(wrap(res_tmp));
  return wrap(res);
  
}

/***R
library(igraph)
library(ergm)
set.seed(1)
m <- rbernoulli(5, p = .3)
w <- geodesic(m)
w[w<0] <- 0
w2 <- geodist(m)
w2$gdist[!is.finite(w2$gdist)] <- 0

w - w2$gdist


net <- as.network(m)
ig  <- graph_from_adjacency_matrix(m)
microbenchmark::microbenchmark(
  ergmito = geodesic(m),
  sna     = geodist(net),
  igraph  = distances(ig),
  unit = "relative"
)
*/