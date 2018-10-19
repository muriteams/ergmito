#include <Rcpp.h>
using namespace Rcpp;

typedef std::vector< int > vecint;

template<typename Tm>
void contingency_matrix(std::vector<int> & table, const Tm & M1, const Tm & M2) {
  
  // Catching error
  if ((M1.ncol() != M2.ncol()) | (M1.nrow() != M2.nrow())) 
    stop("`a` and `b` must have the same dimensions.");
  
  if (table.size() != 4)
    stop("`table` must be of size 4.");
  
  double ans;
  std::fill(table.begin(), table.end(), 0);
  
  int n = M1.nrow();
  int m = M1.ncol();
  
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      if (i == j) continue;
      else {
        if ((M1(i,j) == M2(i,j)) & (M1(i,j) == 0)) table[0]++;
        else if ((M1(i,j) != M2(i,j)) & (M1(i,j) == 0)) table[1]++;
        else if ((M1(i,j) != M2(i,j)) & (M1(i,j) == 1)) table[2]++;
        else if ((M1(i,j) == M2(i,j)) & (M1(i,j) == 1)) table[3]++;
      }
  
}


template<typename Tm> inline
std::vector<int> contingency_matrix(const Tm & M1, const Tm & M2) {
  
  std::vector<int> table(4);
  contingency_matrix< Tm >(table, M1, M2);
  
  return table;
  
}

double starwid(
    const IntegerMatrix & R1,
    const IntegerMatrix & R2,
    bool normalized = false
  ) {
  
  std::vector<int> table = contingency_matrix(R1, R2);
  int n = R1.nrow();
  return (n*table[0] - (table[0] + table[1])*(table[0] + table[2]))/
          (n*table[0] + (table[0] + table[1])*(table[0] + table[2]));
  
}

double S14(
    const IntegerMatrix & R1,
    const IntegerMatrix & R2,
    bool normalized = false
  ) {
  
  std::vector<int> table = contingency_matrix(R1, R2);
  
  return sqrt(
    (
        (table[0] + 1e-15)/(table[0]+table[2]+1e-15) - (table[1] + 1e-15)/(table[1]+table[3]+1e-15))*
          ((table[0] + 1e-15)/(table[0]+table[1]+1e-15) - (table[2] + 1e-15)/(table[2]+table[3]+1e-15))
  );
        
}

double hamming(
    const IntegerMatrix & R1,
    const IntegerMatrix & R2,
    bool normalized = false
  ) {
  
  double ans = 0;
  int n = R1.nrow(), m = R1.ncol();
  for (int i = 0; i<n; ++i) 
    for (int j = 0; j<m; ++j) {
      if (i == j)
        continue;
      
      if (R1(i,j) != R2(i,j))
        ans++;
      
    }
  
  if (normalized) {
    double dn = (double) n;
      ans /= (dn*(dn - 1.0)) ;
  }
    
  return ans;
  
}

double dennis(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    double normalized = false
) {
  
  std::vector<int> table = contingency_matrix(M1, M2);
  
  int n = M1.nrow();
  
  return (table[0]*table[3] - table[1]*table[2])/
    sqrt(n*(table[0] + table[1])*(table[0] + table[2]));
  
}

// -----------------------------------------------------------------------------

NumericMatrix allsimilarities(
    const ListOf<IntegerMatrix> & M,
    bool normalized = false,
    std::function<double(IntegerMatrix, IntegerMatrix, bool)> f = NULL
) {
  
  int N = M.size();
  int NN = N*(N-1)/2;
  NumericVector ans(NN);
  NumericVector I(NN),J(NN);
  
  int pos = 0;
  for (int i = 0; i < N; ++i)
    for (int j = i; j < N; ++j)
      if (i == j) continue;
      else {
        I.at(pos) = i + 1;
        J.at(pos) = j + 1;
        ans.at(pos++) = f(M[i], M[j], normalized);
      }
      
      return cbind(I,J,ans);
  
}

// [[Rcpp::export(name=".S14")]]
NumericMatrix S14_(const ListOf<IntegerMatrix> & M) {
  return allsimilarities(M, false, S14);
  }

// [[Rcpp::export(name=".hamming")]]
NumericMatrix hamming_(const ListOf<IntegerMatrix> & M, bool normalized=false) {
  return allsimilarities(M, false, hamming);
  }

// [[Rcpp::export(name=".starwid")]]
NumericMatrix starwid_(const ListOf<IntegerMatrix> & M, bool normalized=false) {
  return allsimilarities(M, false, starwid);
}

// [[Rcpp::export(name=".dennis")]]
NumericMatrix dennis_(const ListOf<IntegerMatrix> & M, bool normalized=false) {
  return allsimilarities(M, false, dennis);
}
