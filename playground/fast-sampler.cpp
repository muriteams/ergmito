#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]



// [[Rcpp::export]]
arma::vec fast_sampler(int n, const arma::vec & weights) {
  
  arma::vec W = arma::cumsum(weights);
  arma::vec ans(n);
  
  int N       = weights.size();
  double half = N/2;
  int halfi   = floor(half);
  int lb=1, ub=N;
  
  NumericVector S = runif(n);
  W /= W.at(N-1);
  
  int i = 0;
  double dif_left, dif_right;
  // 
  // Rprintf("half: %i\n", half);
  // print(S);
  // W.print();
  
  for (NumericVector::const_iterator it = S.begin(); it != S.end(); ++it) {
    
    int pos = halfi;
    int j   = -1;
    lb=1, ub=N;
    while (j++ < 50) {
      
      // Have we reached the limit?
      if (pos < 1) {
        ans[i++] = 0;
        break;
      } else if (pos > (N - 1)) {
        ans[i++] = N - 1;
        break;
      }
      
      // Checking if is in the range
      dif_left  = *it - W[pos - 1];
      dif_right = W[pos] - *it;
      
      // if (j > 45)
        // Rprintf("W[pos]: %.2f, W[pos - 1]: %.2f, *it: %.2f, pos: %i lb: %i, ub: %i\n", W[pos], W[pos - 1], *it, pos, lb, ub);
      
      // Almost there
      if ((dif_left >= 0.0) & (dif_right >= 0.0)) {
        
        Rprintf("W[pos]: %.2f, W[pos - 1]: %.2f, *it: %.2f, pos: %i lb: %i, ub: %i, j:%i\n", W[pos], W[pos - 1], *it, pos, lb, ub, j);
        
        if (dif_left > dif_right)
          ans[i++] = pos;
        else
          ans[i++] = pos - 1;
        
        break;
      }
      
      // Updating j
      if (*it < W[pos]) {
        ub  = pos;
        pos -= fmax(1, floor((pos - lb)/2));
      } else {
        lb  = pos;
        pos += fmax(1, floor((ub - pos)/2));
      }
      
    }
    
  }
  
  return ans;
  
}
