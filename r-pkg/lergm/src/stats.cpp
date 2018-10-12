#include <Rcpp.h>
using namespace Rcpp;

double s14(const IntegerMatrix & R1, const IntegerMatrix & R2) {
  
  // Catching error
  if ((R1.ncol() != R2.ncol()) | (R1.nrow() != R2.nrow())) 
    stop("`a` and `b` must have the same dimensions.");
  
  double ans;
  double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
  int n = R1.nrow();
  int m = R1.ncol();
  
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      if (i == j) continue;
      else {
        if ((R1(i,j) == R2(i,j)) & (R1(i,j) == 0)) a++;
        else if ((R1(i,j) != R2(i,j)) & (R1(i,j) == 1)) c++;
        else if ((R1(i,j) != R2(i,j)) & (R1(i,j) == 0)) b++;
        else if ((R1(i,j) == R2(i,j)) & (R1(i,j) == 1)) d++;
      }
      
  
  return sqrt(
    ((a + 1e-15)/(a+c+1e-15) - (b + 1e-15)/(b+d+1e-15))*((a + 1e-15)/(a+b+1e-15) - (c + 1e-15)/(c+d+1e-15))
  );
        
}

// [[Rcpp::export(name=".S14")]]
NumericMatrix S14(const std::vector<IntegerMatrix> & M) {
  
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
        ans.at(pos++) = s14(M.at(i), M.at(j));
      }
  
  return cbind(I,J,ans);
  
}

double hamming(const IntegerMatrix & R1, const IntegerMatrix & R2) {
  
  double ans = 0;
  int n = R1.nrow(), m = R1.ncol();
  for (int i = 0; i<n; ++i) 
    for (int j = 0; j<m; ++j) {
      
      if (i == j)
        continue;
      
      if (R1(i,j) == R2(i,j))
        ans++;
      
    }
  
    
  return ans;
  
}

