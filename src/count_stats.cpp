#include <Rcpp.h>
using namespace Rcpp;

/* This function may be useful in the future:
 * List x(1);
 * x.containsElementNamed("Casa")
 */

inline double count_mutual(const IntegerMatrix & x, const NumericVector & A) {
  
  int count = 0;
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = i; j < x.nrow(); ++j)
      if (i != j && x(i,j) + x(j, i) > 1)
        ++count;

#ifdef ERGMITO_COUNT_STATS_DEBUG
  print(x);
  Rprintf("[debug count_mutual] %d\n", count);  
#endif
  
  return (double) count;
}

inline double count_edges(const IntegerMatrix & x, const NumericVector & A) {
  
  int count = 0;
  for (IntegerMatrix::const_iterator it = x.begin(); it != x.end(); ++it)
    if (*it > 0)
      ++count;
    
#ifdef ERGMITO_COUNT_STATS_DEBUG
    print(x);
    Rprintf("[debug count_edges] %d\n", count);
#endif
      
  return (double) count;
}

inline double count_ttriad(const IntegerMatrix & x, const NumericVector & A) {
  
  int count = 0;
  
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = 0; j < x.nrow(); ++j) {
      
      if (x(i,j) == 0)
        continue;
      
      for (int k = 0; k < x.nrow(); ++k) {
        
        // Label 1
        if (x(i,j) == 1 && x(i,k) == 1 && x(j,k) == 1)
          // if (x(j, i) == 0 && x(k,i) == 0 && x(k,j) == 0)
          ++count;
        
      }
    }
  
  return (double) count;
  
}

inline double count_ctriad(const IntegerMatrix & x, const NumericVector & A) {
  
  int count = 0;
  
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = 0; j < i; ++j) {
      
      if (x(i,j) == 0)
        continue;
      
      for (int k = 0; k < i; ++k) {
        
        // Label 1
        if (x(i, j) == 1 && x(j, k) == 1 && x(k, i) == 1)
          // if (x(j, i) == 0 && x(k, j) == 0 && x(i, k) == 0)
            ++count;
          
      }
    }
    
    return (double) count;
  
}

inline double count_nodecov(const IntegerMatrix & x, const NumericVector & A, bool ego) {

  double count = 0.0;
  
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = 0; j < x.nrow(); ++j)
      if (x(i,j) == 1)
        count += A[ego ? i : j];
  
  return count;

}

inline double count_nodeicov(const IntegerMatrix & x, const NumericVector & A) {
  return count_nodecov(x, A, false);
}

inline double count_nodeocov(const IntegerMatrix & x, const NumericVector & A) {
  return count_nodecov(x, A, true);
}


inline double count_nodematch(const IntegerMatrix & x, const NumericVector & A) {
  int count = 0;
  
  for (int i = 0; i < x.nrow(); ++i)
    for (int j = 0; j < x.nrow(); ++j)
      if (x(i,j) == 1 && A.at(i) == A.at(j))
        ++count;
      
      return (double) count;
}

inline double count_triangle(const IntegerMatrix & x, const NumericVector & A) {
  
  return count_ctriad(x, A) + count_ttriad(x, A);
  
}

typedef double (*ergm_term_fun)(const IntegerMatrix & x, const NumericVector & A);

void get_ergm_term(std::string term, ergm_term_fun & fun) {
  
  if (term == "mutual")         fun = &count_mutual;
  else if (term == "edges")     fun = &count_edges;
  else if (term == "ttriad")    fun = &count_ttriad;
  else if (term == "ctriad")    fun = &count_ctriad;
  else if (term == "nodeicov")  fun = &count_nodeicov;
  else if (term == "nodeocov")  fun = &count_nodeocov;
  else if (term == "nodematch") fun = &count_nodematch;
  else if (term == "triangle")  fun = &count_triangle;
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
  
  bool uses_attributes = false;
  NumericVector A_empty(0);
  if (A[0].size() != 0) {
    if (A.size() != n)
      stop("The number of attributes in `A` differs from the number of adjacency matrices.");
    
    uses_attributes = true;
  }
  
  NumericMatrix ans(n, k);
  ergm_term_fun fun;
  
  for (int j = 0; j < k; ++j) {

    // Getting the function
    get_ergm_term(terms[j], fun);

    if (uses_attributes) {
      for (int i = 0; i < n; ++i)
        ans.at(i, j) = fun(X[i], A[i]);
    } else {
      for (int i = 0; i < n; ++i)
        ans.at(i, j) = fun(X[i], A_empty);
    }
    
  }
//   
// #ifdef ERGMITO_COUNT_STATS_DEBUG
//   print(ans);
// #endif
//   
  return ans;
  
}

/***R

library(ergmito)
x <- powerset(5)

# Problems with large lists, ... again we'll need to do this by chunks.
s <- count_stats(x[1:1e5], c("mutual", "edges"))

*/