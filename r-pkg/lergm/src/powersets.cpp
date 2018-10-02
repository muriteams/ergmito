#include "lergm_types.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
vecint make_sets(int n) {
  
  int m = n*(n-1);
  vecint ans(m);
  
  int pos=-1, k=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      if (i!=j)
        ans[++pos] = k;
      ++k;
    }
    
    return ans;
}

// [[Rcpp::export(name=.powerset)]]
vecvecint powerset(int n, bool force = false) {
  
  if (n > 5 && !force)
    Rcpp::stop("In order to generate power sets for n>5 force must be set to `TRUE`.");
  
  int m = n*(n-1);
  vecint set = make_sets(n);
  vecvecint sets(pow(2.0, m));
  vecint v(1);
  
  int j = 0,k;
  for (int i=0; i<m; i++) {
    k = j;
    for (int s = 0; s < k; s++) {
      sets.at(j) = sets.at(s);
      sets.at(j++).push_back(set[i]);
    }
    v.at(0) = set[i];
    sets.at(j++) = v;
  }

  return sets;
  
}


