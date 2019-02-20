#include <Rcpp.h>
using namespace Rcpp;

inline double count_mutual(const IntegerMatrix & x, const NumericVector & A) {
  
  double count = 0.0;
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = i; j < x.nrow(); ++j)
      if (i != j && x(i,j) + x(j, i) > 1)
        count += 1.0;
  
  return count;
}

inline double count_edges(const IntegerMatrix & x, const NumericVector & A) {
  
  double count = 0.0;
  
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = 0; j < x.nrow(); ++j)
      if (x(i,j) > 0)
        count += 1;
      
  return count;
}

inline double count_ttriad(const IntegerMatrix & x, const NumericVector & A) {
  
  double count = 0.0;
  
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = 0; j < x.nrow(); ++j) {
      
      if (x(i,j) == 0)
        continue;
      
      for (int k = 0; k < x.nrow(); ++k) {
        
        // Label 1
        if (x(i,j) == 1 && x(i,k) == 1 && x(j,k) == 1)
          // if (x(j, i) == 0 && x(k,i) == 0 && x(k,j) == 0)
          count += 1.0;
        
      }
    }
  
  return count;
  
}

inline double count_ctriad(const IntegerMatrix & x, const NumericVector & A) {
  
  double count = 0.0;
  
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = 0; j < i; ++j) {
      
      if (x(i,j) == 0)
        continue;
      
      for (int k = 0; k < i; ++k) {
        
        // Label 1
        if (x(i, j) == 1 && x(j, k) == 1 && x(k, i) == 1)
          // if (x(j, i) == 0 && x(k, j) == 0 && x(i, k) == 0)
            count += 1.0;
          
      }
    }
    
    return count;
  
}

inline double count_nodecov(const IntegerMatrix & x, const NumericVector & A, bool ego) {

  double count = 0.0;
  
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = 0; j < x.nrow(); ++j)
      if (x(i,j) == 1)
        count += A.at(ego ? i : j);
  
  return count;

}

inline double count_nodeicov(const IntegerMatrix & x, const NumericVector & A) {
  return count_nodecov(x, A, false);
}

inline double count_nodeocov(const IntegerMatrix & x, const NumericVector & A) {
  return count_nodecov(x, A, true);
}


inline double count_nodematch(const IntegerMatrix & x, const NumericVector & A) {
  double count = 0.0;
  
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = 0; j < x.nrow(); ++j)
      if (x(i,j) == 1 && A.at(i) == A.at(j))
        count += 1.0;
      
      return count;
}

inline double count_triangle(const IntegerMatrix & x, const NumericVector & A) {
  
  return count_ctriad(x, A) + count_ttriad(x, A);
  
}

typedef double (*ergm_term_fun)(const IntegerMatrix & x, const NumericVector & A);

void get_ergm_term(std::string term, ergm_term_fun & fun) {
  
  if (term == "mutual")      fun = &count_mutual;
  else if (term == "edges")  fun = &count_edges;
  else if (term == "ttriad") fun = &count_ttriad;
  else if (term == "ctriad") fun = &count_ctriad;
  else if (term == "nodeicov") fun = &count_nodeicov;
  else if (term == "nodeocov") fun = &count_nodeocov;
  else if (term == "nodematch") fun = &count_nodematch;
  else if (term == "triangle") fun = &count_triangle;
  else
    stop("The term %s is not available.", term);
  
  return;
  
}

// [[Rcpp::export]]
CharacterVector count_available(int i = 0) {
  return CharacterVector::create(
    "mutual",
    "edges",
    "ttriad",
    "ctriad",
    "nodeicov",
    "nodeocov",
    "nodematch",
    "triangle"
  );
}

// Count Network Statistics
// [[Rcpp::export(name="count_stats.")]]
NumericMatrix count_stats(
    const ListOf< IntegerMatrix > & X,
    const std::vector< std::string > & terms,
    const ListOf< NumericVector > & A
  ) {
  
  
  int n = X.size();
  int k = terms.size();
  
  NumericMatrix ans(n, k);
  ergm_term_fun fun;
  
  int i = 0;
  for (int j = 0; j < k; ++j) {

    // Getting the function
    get_ergm_term(terms.at(j), fun);

    for (ListOf< IntegerMatrix >::const_iterator x = X.begin(); x != X.end(); ++x)
      ans(i++, j) = fun(*x, A[x.index()]);
    
    i = 0;
    
  }
  
  return ans;
  
}

/***R

library(ergmito)
x <- powerset(5)

# Problems with large lists, ... again we'll need to do this by chunks.
s <- count_stats(x[1:1e5], c("mutual", "edges"))

*/