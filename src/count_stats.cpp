#include <RcppArmadillo.h>
using namespace Rcpp;

/* This function may be useful in the future:
 * List x(1);
 * x.containsElementNamed("Casa")
 */

inline double count_mutual(const IntegerMatrix & x, const NumericVector & A) {
  
  unsigned int count = 0u;
  for (unsigned int i = 0u; i < x.nrow(); ++i)
    for (unsigned int j = i; j < x.nrow(); ++j)
      if (i != j && x(i,j) + x(j, i) > 1)
        ++count;

#ifdef ERGMITO_COUNT_STATS_DEBUG
  print(x);
  Rprintf("[debug count_mutual] %d\n", count);  
#endif
  
  return (double) count;
}

inline double count_edges(const IntegerMatrix & x, const NumericVector & A) {
  
  unsigned int count = 0u;
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
  
  unsigned int count = 0u;
  
  for (unsigned int i = 0u; i < x.nrow(); ++i)
    for (unsigned int j = 0u; j < x.nrow(); ++j) {
      
      if (x(i,j) == 0)
        continue;
      
      for (unsigned int k = 0u; k < x.nrow(); ++k) {
        
        // Label 1
        if (x(i,j) == 1 && x(i,k) == 1 && x(j,k) == 1)
          ++count;
        
      }
    }
  
  return (double) count;
  
}

inline double count_ctriad(const IntegerMatrix & x, const NumericVector & A) {
  
  unsigned int count = 0u;
  
  for (unsigned int i = 0u; i < x.nrow(); ++i)
    for (unsigned int j = 0u; j < i; ++j) {
      
      if (x(i,j) == 0)
        continue;
      
      for (unsigned int k = 0u; k < i; ++k) {
        
        // Label 1
        if (x(i, j) == 1 && x(j, k) == 1 && x(k, i) == 1)
          ++count;
          
      }
    }
    
    return (double) count;
  
}

inline double count_nodecov(const IntegerMatrix & x, const NumericVector & A, bool ego) {

  double count = 0.0;
  
  for (unsigned int i = 0u; i < x.nrow(); ++i)
    for (unsigned int j = 0u; j < x.nrow(); ++j)
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
  unsigned int count = 0u;
  
  for (unsigned int i = 0u; i < x.nrow(); ++i)
    for (unsigned int j = 0u; j < x.nrow(); ++j)
      if (x(i,j) == 1 && A.at(i) == A.at(j))
        ++count;
      
      return (double) count;
}

inline double count_triangle(const IntegerMatrix & x, const NumericVector & A) {
  
  return count_ctriad(x, A) + count_ttriad(x, A);
  
}

inline double count_balance(const IntegerMatrix & x, const NumericVector & A) {
  
  unsigned int count = 0u, n = x.nrow();
  int s;
  
  for (unsigned int i = 0u; i < n; ++i)
    for (unsigned int j = 0u; j < n; ++j) {
  
      if (i == j)
        continue;
  
      // We compute this statistic for selecting two possible cases: disconnected
      // and mutual ties. The two relevant cases for balanced triads.
      s = x(i, j) + x(j, i);

      if (s == 1) 
        continue;
      
      // Case in which i and j are disconnected. If these two are disconnected,
      // then the loop through k should be truncated as a function of i since
      // otherwise we will be double counting.
      else if (s == 0) { // Triad 102
        
        for (unsigned int k = 0u; k < i; ++k) {
          
          if (k == j)
            continue;
        
          if (x(i, k) == 0 || x(k, i) == 0 || x(j, k) == 1 || x(k, j) == 1)
            continue;
          
          ++count;
        

        }
        
        
        // Case in which they have a mutual tie, then we also have to truncate the
        // loop over j, otherwise triple counting.
      } else if (s == 2) { // Triad 300
        
        if (j >= i)
          break;
        
        for (unsigned int k = 0u; k < j; ++k) {
          
          if (x(i, k) == 0 || x(k, i) == 0 || x(j, k) == 0 || x(k, j) == 0)
            continue;
          
          ++count;
        }
        
      }
      
    }
  
  return (double) count;
  
}

// Triadic census --------------------------------------------------------------
inline double count_t300(const IntegerMatrix & x, const NumericVector & A) {

  unsigned int count = 0u, n = x.nrow();

  for (unsigned int i = 0u; i < n; ++i) {
    
    for (unsigned int j = 0u; j < i; ++j) {
      
      if (x(i, j) == 0 || x(j, i) == 0)
        continue;
      
      for (unsigned int k = 0u; k < j; ++k) {
        
        if (x(i, k) == 0 || x(k, i) == 0 || x(j, k) == 0 || x(k, j) == 0)
          continue;
        
        ++count;
        
      }
      
    }
    
  }
  
  return (double) count;
}

inline double count_t102(const IntegerMatrix & x, const NumericVector & A) {
  unsigned int count = 0u, n = x.nrow();
  
  for (unsigned int i = 0u; i < n; ++i) {
    
    for (unsigned int j = 0u; j < n; ++j) {
      
      if (x(i, j) == 1 || x(j, i) == 1)
        continue;
      
      for (unsigned int k = 0u; k < i; ++k) {
        
        if (x(i, k) == 0 || x(k, i) == 0 || x(j, k) == 1 || x(k, j) == 1)
          continue;
        
        ++count;
        
      }
      
    }
    
  }
  
  return (double) count;
}

// Callers ---------------------------------------------------------------------
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
  else if (term == "balance")   fun = &count_balance;
  else if (term == "t300")      fun = &count_t300;
  else if (term == "t102")      fun = &count_t102;
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
    "triangle",
    "balance",
    "t300",
    "t102"
  );
}

// Count Network Statistics
// [[Rcpp::export(name="count_stats.")]]
NumericMatrix count_stats(
    const ListOf< IntegerMatrix > & X,
    const std::vector< std::string > & terms,
    const ListOf< NumericVector > & A
  ) {
  
  // List b;
  // b.containsElementNamed("parameter1");
  
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
  
  for (unsigned int j = 0u; j < k; ++j) {
      
    // Getting the function
    get_ergm_term(terms[j], fun);

    if (uses_attributes) {
      for (unsigned int i = 0u; i < n; ++i) {
        
        // Checking dimmensions
        if (X[i].nrow() != X[i].ncol())
          stop("Matrix %i is not a square matrix.", i);
        
        ans.at(i, j) = fun(X[i], A[i]);
      }
        
    } else {
      for (unsigned int i = 0u; i < n; ++i) {
        // Checking dimmensions
        if (X[i].nrow() != X[i].ncol())
          stop("Matrix %i is not a square matrix.", i);
        
        ans.at(i, j) = fun(X[i], A_empty);
      }
        
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

inline IntegerMatrix geodesici(const arma::imat & x, bool force = false) {
  
  int n = x.n_rows;
  if (n != (int) x.n_cols)
    stop("Not a square matrix.");
  
  if (n > 100 && !force)
    stop("This is not the best way for computing distances for n > 100 (see ?geodesic).");
  
  arma::imat res(x.n_rows, x.n_cols);
  res.fill(-1);
  arma::imat res_tmp = x;
  
  // List of indices to check
  unsigned int changes = 0u;
  for (unsigned int iter = 0u; iter < n; ++iter) {
    changes = 0;
    for (unsigned int i = 0u; i < n; ++i)
      for (unsigned int j = 0u; j < n; ++j) {
        if (i != j && res_tmp.at(i, j) != 0 && (res.at(i, j) == -1)) 
          res.at(i, j) = iter + 1, ++changes;
      }
      
      // Powermatrix
      res_tmp *= x;
  }
  
  return wrap(res);
  
}

// [[Rcpp::export(name = "geodesic.")]]
ListOf< IntegerMatrix > geodesic(
    const std::vector< IntegerMatrix > & X,
    bool force = false) {
  
  ListOf< IntegerMatrix > res(X.size());
  
  unsigned int nX = X.size();
  for (unsigned int i = 0u; i < nX; ++i)
    res[i] = geodesici(as< arma::imat >(X[i]), force);
  
  return res;
  
  
}
