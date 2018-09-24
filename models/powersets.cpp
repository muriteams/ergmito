#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector< int > make_sets(int n) {
  
  int m = n*(n-1);
  std::vector< int > ans(m);
  
  int pos=-1, k=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      if (i!=j)
        ans[++pos] = k;
      ++k;
    }
    
    return ans;
}

// [[Rcpp::export]]
std::vector< std::vector<int> > powerset(int n) {
  
  int m = n*(n-1);
  std::vector< int > set = make_sets(n);
  std::vector< std::vector< int > > sets(pow(2.0, m));
  std::vector< int > v(1);
  
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


