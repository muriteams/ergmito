#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector< std::vector<int> > powerset(int n) {
  
  int m = n*(n-1);
  std::vector< std::vector< int > > sets(pow(2, m));
  std::vector< int > v(1);
  
  int j = 0,k;
  for (int i=0; i<m; i++) {
    k = j;
    for (int s = 0; s < k; s++) {
      sets.at(j) = sets.at(s);
      sets.at(j++).push_back(i);
    }
    v.at(0) = i;
    sets.at(j++) = v;
  }

  return sets;
  
}