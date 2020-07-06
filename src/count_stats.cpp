#include <RcppArmadillo.h>
using namespace Rcpp;

/* This function may be useful in the future:
 * List x(1);
 * x.containsElementNamed("Casa")
 */

inline double count_mutual(const IntegerMatrix & x, const NumericVector & A) {
  
  unsigned int count = 0u, n = (unsigned int) x.nrow();
  for (unsigned int i = 0u; i < n; ++i)
    for (unsigned int j = i; j < n; ++j)
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
  unsigned int n = (unsigned int) x.nrow();
  
  for (unsigned int i = 0u; i < n; ++i)
    for (unsigned int j = 0u; j < n; ++j) {
      
      if (x(i,j) == 0)
        continue;
      
      for (unsigned int k = 0u; k < n; ++k) {
        
        // Label 1
        if (x(i,j) == 1 && x(i,k) == 1 && x(j,k) == 1)
          ++count;
        
      }
    }
  
  return (double) count;
  
}

inline double count_ctriad(const IntegerMatrix & x, const NumericVector & A) {
  
  unsigned int count = 0u;
  unsigned int n = (unsigned int) x.nrow();
  
  for (unsigned int i = 0u; i < n; ++i)
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

inline double count_absdiff(const IntegerMatrix & x, const NumericVector & A) {
  double count = 0.0;
  unsigned int n = (unsigned int) x.nrow();
  for (unsigned int i = 0u; i < n; ++i)
    for (unsigned int j = 0u; j < n; ++j)
      if (x(i,j) == 1)
        count += fabs(A[i] - A[j]);
  
  return count;
}

inline double count_nodecov(const IntegerMatrix & x, const NumericVector & A, bool ego) {

  double count = 0.0;
  unsigned int n = (unsigned int) x.nrow();
  
  for (unsigned int i = 0u; i < n; ++i)
    for (unsigned int j = 0u; j < n; ++j)
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
  
  unsigned int count = 0u, n = (unsigned int) x.nrow();
  
  for (unsigned int i = 0u; i < n; ++i)
    for (unsigned int j = 0u; j < n; ++j)
      if (x(i,j) == 1 && A.at(i) == A.at(j))
        ++count;
      
  return (double) count;
  
}

inline double count_triangle(const IntegerMatrix & x, const NumericVector & A) {
  
  return count_ctriad(x, A) + count_ttriad(x, A);
  
}

inline double count_idegree15(const IntegerMatrix & x, const NumericVector & A) {
  return sum(pow(colSums(x), 1.5));
}
inline double count_odegree15(const IntegerMatrix & x, const NumericVector & A) {
  return sum(pow(rowSums(x), 1.5));
}
// inline double count_degree15(const IntegerMatrix & x, const NumericVector & A) {
//   return sum(pow(rowSums(x) + colSums(x), 1.5));
// }

inline double count_star1(
    const IntegerMatrix & x,
    const NumericVector & A,
    bool out = true
) {
  
  // We only use the first element of the count
  int count = 0;
  unsigned int i, j, n = (unsigned int) x.nrow();
  
  // If it is outgoing, then point to it...
  unsigned int * head_j,  * tail_j;
  if (out) 
    head_j = &i, tail_j = &j;
  else 
    head_j = &j, tail_j = &i;
  
  for (i = 0u; i < n; ++i) 
    for (j = 0u; j < n; ++j) {
      
      if (i == j)
        continue;
      
      if ((x(*head_j, *tail_j) == 1) &&
          (A.size() == 0 || (A[i] == A[j])))
        count += 1;
    }
    
  return (double) count;
  
}

inline double count_istar1(const IntegerMatrix & x, const NumericVector & A) {
  
  return count_star1(x, A, false);
  
}

inline double count_ostar1(const IntegerMatrix & x, const NumericVector & A) {
  
  return count_star1(x, A, true);
  
}

inline double count_star2(
    const IntegerMatrix & x,
    const NumericVector & A,
    bool out = true
  ) {
  
  // We only use the first element of the count
  int count = 0;
  unsigned int i, j, k, n = (unsigned int) x.nrow();
  
  // If it is outgoing, then point to it...
  unsigned int * head_j,  * tail_j, * head_k, * tail_k;
  if (out) 
    head_j = &i, tail_j = &j, head_k = &i, tail_k = &k;
  else 
    head_j = &j, tail_j = &i, head_k = &k, tail_k = &i;
  
  for (i = 0u; i < n; ++i) 
    for (j = 0u; j < n; ++j) {
    
      if (i == j)
        continue;
      
      for (k = j; k < n; ++k) {
        
        if ((i == k) | (k == j))
          continue;
        
        if ((x(*head_j, *tail_j) == 1) && (x(*head_k, *tail_k) == 1) &&
            (A.size() == 0 || (A[i] == A[j] && A[i] == A[k])))
          count += 1;
      }
    }
      
  return (double) count;
  
}

inline double count_istar2(const IntegerMatrix & x, const NumericVector & A) {
  
  return count_star2(x, A, false);
  
}

inline double count_ostar2(const IntegerMatrix & x, const NumericVector & A) {
  
  return count_star2(x, A, true);
  
}

inline double count_star3(const IntegerMatrix & x, const NumericVector & A,
                          bool out) {
  
  // We only use the first element of the count
  int count = 0;
  unsigned int i, j, k, l, n = (unsigned int) x.nrow();
  
  // If it is outgoing, then point to it...
  unsigned int * head_j,  * tail_j, * head_k, * tail_k, * head_l, * tail_l;
  if (out) 
    head_j = &i, tail_j = &j,
      head_k = &i, tail_k = &k,
      head_l = &i, tail_l = &l;
  else 
    head_j = &j, tail_j = &i,
      head_k = &k, tail_k = &i,
      head_l = &l, tail_l = &i;
  
  for (i = 0u; i < n; ++i)
    for (j = 0u; j < n; ++j) {
      
      if (i == j)
        continue;
      
      for (k = j; k < n; ++k) {
        
        if ((i == k) | (j == k))
          continue;
        
        for (l = k; l < n; ++l) {
          
          if ((i == l) | (k == l) | (j == l))
            continue;
          
          if (
              (x(*head_j, *tail_j) == 1) && (x(*head_k, *tail_k) == 1) && (x(*head_l, *tail_l) == 1) &&
                ((A.size() == 0) || (A[i] == A[j] && A[i] == A[k] && A[i] == A[l]))
          )
            count += 1;
        }
          
      }
    }
        
  return (double) count;
        
}

inline double count_ostar3(const IntegerMatrix & x, const NumericVector & A) {
  return count_star3(x, A, true);
}

inline double count_istar3(const IntegerMatrix & x, const NumericVector & A) {
  return count_star3(x, A, false);
}

inline double count_star4(const IntegerMatrix & x, const NumericVector & A,
                           bool out) {
  
  // We only use the first element of the count
  int count = 0;
  unsigned int i, j, k, l, m, n = (unsigned int) x.nrow();
  
  // If it is outgoing, then point to it...
  unsigned int * head_j,  * tail_j, * head_k, * tail_k, * head_l, * tail_l,
  *head_m, *tail_m;
  if (out) 
    head_j = &i, tail_j = &j,
      head_k = &i, tail_k = &k,
      head_l = &i, tail_l = &l,
      head_m = &i, tail_m = &m;
  else 
    head_j = &j, tail_j = &i,
      head_k = &k, tail_k = &i,
      head_l = &l, tail_l = &i,
      head_m = &m, tail_m = &i;
  
  for (i = 0u; i < n; ++i) 
    for (j = 0u; j < n; ++j) { 
      
      if (i == j)
        continue;
      
      for (k = j; k < n; ++k) {
        
        if ((i == k) | (j == k))
          continue;
        
        for (l = k; l < n; ++l) {
          
          if ((i == l) | (j == l) | (k == l))
            continue;
          
          for (m = l; m < n; ++m) {
            
            if ((i == m) | (j == m) | (k == m) | (l == m))
              continue;
            
            if (
                (x(*head_j, *tail_j) == 1) && (x(*head_k, *tail_k) == 1) && (x(*head_l, *tail_l) == 1) && (x(*head_m, *tail_m) == 1) &&
                  ((A.size() == 0u) || (A[i] == A[j] && A[i] == A[k] && A[i] == A[l] && A[i] == A[m]))
              )
              count += 1;
          }
        }
      }
    }
          
  return (double) count;
          
}

inline double count_ostar4(const IntegerMatrix & x, const NumericVector & A) {
  return count_star4(x, A, true);
}
  
inline double count_istar4(const IntegerMatrix & x, const NumericVector & A) {
  return count_star4(x, A, false);
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
  
  if (term == "mutual")          fun = &count_mutual;
  else if (term == "edges")      fun = &count_edges;
  else if (term == "ttriad")     fun = &count_ttriad;
  else if (term == "ctriad")     fun = &count_ctriad;
  else if (term == "ctriple")    fun = &count_ctriad;
  else if (term == "nodeicov")   fun = &count_nodeicov;
  else if (term == "nodeocov")   fun = &count_nodeocov;
  else if (term == "nodematch")  fun = &count_nodematch;
  else if (term == "triangle")   fun = &count_triangle;
  else if (term == "balance")    fun = &count_balance;
  else if (term == "t300")       fun = &count_t300;
  else if (term == "t102")       fun = &count_t102;
  else if (term == "absdiff")    fun = &count_absdiff;
  // else if (term == "degree1.5")  fun = &count_degree15;
  else if (term == "idegree1.5") fun = &count_idegree15;
  else if (term == "odegree1.5") fun = &count_odegree15;
  else if (term == "ostar1")     fun = &count_ostar1;
  else if (term == "ostar2")     fun = &count_ostar2;
  else if (term == "ostar3")     fun = &count_ostar3;
  else if (term == "ostar4")     fun = &count_ostar4;
  else if (term == "istar1")     fun = &count_istar1;
  else if (term == "istar2")     fun = &count_istar2;
  else if (term == "istar3")     fun = &count_istar3;
  else if (term == "istar4")     fun = &count_istar4;
  else
    stop("The term %s is not available in ergmito.", term);
  
  return;
  
}

// [[Rcpp::export(rng = false)]]
std::vector< std::string > count_available(int i = 0) {
  return {
    "mutual",
    "edges",
    "ttriad",
    "ctriad", "ctriple",
    "nodeicov",
    "nodeocov",
    "nodematch",
    "triangle",
    "balance",
    "t300",
    "t102",
    "absdiff",
    // "degree1.5",
    "idegree1.5",
    "odegree1.5",
    "ostar1", "ostar2", "ostar3", "ostar4",
    "istar1", "istar2", "istar3", "istar4"}
  ;
}

// Count Network Statistics
// [[Rcpp::export(name="count_stats.", rng = false)]]
NumericMatrix count_stats(
    const ListOf< IntegerMatrix > & X,
    const std::vector< std::string > & terms,
    const ListOf< NumericVector > & A
  ) {
  
  // List b;
  // b.containsElementNamed("parameter1");
  
  typedef unsigned int uint;
  uint n = X.size();
  uint k = terms.size();
  
  bool uses_attributes = false;
  bool all_same_attr   = false;
  NumericVector A_empty(0u);
  
  // In this case we are passing the set of attributes to all the data
  // so we need to check that this will be have has expected
  
  uint A0_size = (uint) A[0].size();
  if (A0_size != 0u) {
    
    uint A_size = (uint) A.size();
    
    if (n > 1u) {
      
      // Assuming all have the same attribute
      if (A_size == 1u) {
        
        all_same_attr   = true;
        
        // Need to check that the size fits
        for (auto i = 0u; i < n; ++i)
          if (X[i].nrow() != A[0].size())
            stop(
              "The length of the attributes (%i) does not match the size of one of the matrices (%i).",
              (int) A[0].size(), (int) X[i].nrow()
            );
        
      } else if (A_size != n) 
        stop("The number of attributes in `A` differs from the number of adjacency matrices.");

    } else 
      if (A0_size != (uint) X[0].nrow())
        stop("The length of the attributes does not match the number of vertices.");
      
    uses_attributes = true;
  }
  
  NumericMatrix ans(n, k);
  ergm_term_fun fun;
  
  for (unsigned int j = 0u; j < k; ++j) {
      
    // Getting the function
    get_ergm_term(terms[j], fun);

    if (uses_attributes) {
      for (unsigned int i = 0u; i < n; ++i) {
        
        // Checking dimensions
        if (X[i].nrow() != X[i].ncol())
          stop("Matrix %i is not a square matrix.", i + 1);
        
        if (all_same_attr)
          ans.at(i, j) = fun(X[i], A[0u]);
        else
          ans.at(i, j) = fun(X[i], A[i]);
      }
        
    } else {
      for (unsigned int i = 0u; i < n; ++i) {
        // Checking dimensions
        if (X[i].nrow() != X[i].ncol())
          stop("Matrix %i is not a square matrix.", i + 1);
        
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

inline void geodesici(const arma::imat & x, IntegerMatrix & res, bool force = false) {
  
  unsigned int n = x.n_rows;
  if (n != x.n_cols)
    stop("Not a square matrix.");
  
  if (n > 100u && !force)
    stop("This is not the best way for computing distances for n > 100 (see ?geodesic).");
  
  arma::imat res_tmp = x;
  
  // List of indices to check
  unsigned int nmax = n * 2;
  for (unsigned int iter = 0u; iter < nmax; ++iter) {
    for (unsigned int i = 0u; i < n; ++i) {
      for (unsigned int j = 0u; j < n; ++j) {
        
        if (i == j) {
          res(i, j) = 0;
          continue;
        }
          
        if (i != j && res_tmp.at(i, j) != 0 && (res(i, j) == NA_INTEGER)) 
          res(i, j) = iter + 1;
      }
    }
      
    // Powermatrix
    res_tmp *= x;
  }
  
  return;
  
}

// [[Rcpp::export(name = "geodesic.", rng = false)]]
std::vector< IntegerMatrix > geodesic(
    const std::vector< arma::imat > & X,
    bool force = false
  ) {
  
  std::vector< IntegerMatrix > res;
  res.reserve(X.size());
  
  unsigned int nX = X.size();
  for (unsigned int i = 0u; i < nX; ++i) {
    
    IntegerMatrix tmp(X[i].n_rows, X[i].n_cols);
    tmp.fill(NA_INTEGER);
    res.push_back(tmp);
    geodesici(X[i], res[i], force);
    
  }
  
  return res;
  
  
}
