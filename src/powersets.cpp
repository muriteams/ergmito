#include <numeric>
#include "ergmito_types.h"
#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
vecint make_sets(int n, bool directed = true) {
  
  int N = directed ? n * (n - 1): n * (n - 1) / 2;
  vecint ans(N);
  
  int i, zero = 0;
  int *m = directed ? &zero: &i;
  
  int pos = -1, k = 0;
  for (i = 0; i < n; ++i) {

    for (int j = (*m); j < n; ++j) {
      
      // At the beginning we have to add some to correct for the diagonal
      if (!directed && (i == j))
        k += i;
      
      if (i != j)
        ans[ ++pos ] = k;
      
      ++k;
      
    }
    
  }
    
  return ans;
}

/*I've detected the following problem: When the numner of IntegerVector to be
 * created overpasses a couple of hundred thousand, the stack memory gets filled
 * (overflowed) since by default Rcpp works on the stack instead of the heap.
 * 
 * I have not found yet another solution to this problem, but to divide the
 * creation of IntegerVector into chunks of 200,000.
 */

// This function creates power sets
void powerset(vecvecint * sets, int n, bool directed = true) {
  
  vecint set = make_sets(n, directed);
  int nsets   = set.size();

  int j = 0,k;
  for (int i=0; i < nsets; ++i) {
    k = j;
    for (int s = 0; s < k; ++s) {
      sets->operator[](j) = sets->operator[](s);
      sets->operator[](j++).push_back(set[i]);
    }
    sets->operator[](j++).push_back(set[i]);
  }

  return;
  
}

// [[Rcpp::export(name=".powerset", rng = false)]]
SEXP powerset(int n, bool force = false, bool directed = true) {
  
  if (n > 5 && !force)
    Rcpp::stop("In order to generate power sets for n>5 force must be set to `TRUE`.");
  
  int m = directed? n * (n - 1): n * (n - 1) / 2;
  vecvecint * sets = new vecvecint(pow(2.0, m));

  powerset(sets, n, directed);
  
  return Rcpp::XPtr< vecvecint >(sets, true);
  
}

// // [[Rcpp::export]]
// int print_powerset(SEXP sets) {
//   
//   Rcpp::XPtr< vecvecint > sets_ptr(sets);
//   for (vecvecint::const_iterator it = sets_ptr->begin(); it != sets_ptr->end(); ++it) {
//     print(wrap(*it));
//   }
//   
//   return 0;
//   
// }

// This is another wrapper, this takes care of turning those integer vectors
// into NumericMatrix of size 2.
// [[Rcpp::export(rng = false)]]
List wrap_powerset(
    SEXP sets,
    int from,
    int to,
    int n
) {
  
  // Getting the pointer
  Rcpp::XPtr< vecvecint > sets_ptr(sets);
  
  // Checking boundaries
  int N = sets_ptr->size();
  
  if (from < 0)
    stop("The `from` parameter must be a positive integer.");
  if (to > N)
    stop("The `to` parameter must be smaller than `N`.");
  if (from >= to)
    stop("`from` should be smaller than `to`.");
  
  // Creating a list of objects
  List ans(to - from);
  IntegerVector dim = IntegerVector::create(n, n);
  
  int counter = 0;
  typedef std::vector<int>::const_iterator vintiter;
  // Creating the empty vector
  IntegerVector tmp(n*n);
  // std::cout << from << "\n";
  for (int i = from; i < to; i++) {
    
    tmp.fill(0);
    
    // Filling it with data
    for (vintiter it2 = sets_ptr->operator[](i).begin(); it2 != sets_ptr->operator[](i).end(); ++it2) 
      tmp.at(*it2) = 1;
    
    // Adding the vector to the list, but before set attributes
    tmp.attr("dim") = dim;
    ans.at(counter++) = clone(tmp);
    
  }
  
  return ans;
  
}
