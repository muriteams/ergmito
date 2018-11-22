#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int count_mutual(const IntegerMatrix & x) {
  
  int count = 0;
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = i; j < x.nrow(); ++j)
      if (i != j && x(i,j) + x(j, i) > 1)
        ++count;
  
  return count;
}

// [[Rcpp::export]]
int count_edges(const IntegerMatrix & x) {
  
  int count = 0;
  
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = 0; j < x.nrow(); ++j)
      if (x(i,j) > 0)
        ++count;
      
  return count;
}

typedef int (*ergm_term_fun)(const IntegerMatrix & x);

void get_ergm_term(std::string term, ergm_term_fun & fun) {
  
  if     (term == "mutual") fun = &count_mutual;
  else if (term == "edges") fun = &count_edges;
  else stop("The term %s is not available.", term);
  
  return;
  
}

// [[Rcpp::export]]
IntegerMatrix count_stats(
    const ListOf< IntegerMatrix > & X,
    const std::vector< std::string > & terms
  ) {
  
  
  
  int n = X.size();
  int k = terms.size();
  
  IntegerMatrix ans(n, k);
  ergm_term_fun fun;
  
  int i = 0;
  for (int j = 0; j < k; ++j) {

    // Getting the function
    get_ergm_term(terms.at(j), fun);

    for (ListOf< IntegerMatrix >::const_iterator x = X.begin(); x != X.end(); ++x)
      ans(i++, j) = fun(*x);
    
  }
  
  return ans;
  
}

/***R

library(lergm)
x <- powerset(5)

# Problems with large lists, ... again we'll need to do this by chunks.
s <- count_stats(x[1:1e5], "mutual")

*/