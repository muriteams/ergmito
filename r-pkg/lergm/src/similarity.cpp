#include <Rcpp.h>
using namespace Rcpp;

#define a 0
#define b 1
#define c 2
#define d 3

typedef std::vector< int > vecint;

template<typename Ti, typename Tm>
void contingency_matrix(std::vector<Ti> & table, const Tm & M1, const Tm & M2) {
  
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
        if ((M1(i,j) == M2(i,j)) & (M1(i,j) == 0)) table[a]++;
        else if ((M1(i,j) != M2(i,j)) & (M1(i,j) == 0)) table[b]++;
        else if ((M1(i,j) != M2(i,j)) & (M1(i,j) == 1)) table[c]++;
        else if ((M1(i,j) == M2(i,j)) & (M1(i,j) == 1)) table[d]++;
      }
  
}


template<typename Ti, typename Tm> inline
std::vector<Ti> contingency_matrix(const Tm & M1, const Tm & M2) {
  
  std::vector<Ti> table(4);
  contingency_matrix< Ti, Tm >(table, M1, M2);
  
  return table;
  
}

double starwid(
    const IntegerMatrix & R1,
    const IntegerMatrix & R2,
    bool normalized = false
  ) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(R1, R2);
  int n = R1.nrow();
  return (n*table[a] - (table[a] + table[b])*(table[a] + table[c]))/
          (n*table[a] + (table[a] + table[b])*(table[a] + table[c]));
  
}

double S14(
    const IntegerMatrix & R1,
    const IntegerMatrix & R2,
    bool normalized = false
  ) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(R1, R2);
  
  return sqrt(
    (
        (table[a] + 1e-15)/(table[a]+table[c]+1e-15) - (table[b] + 1e-15)/(table[b]+table[d]+1e-15))*
          ((table[a] + 1e-15)/(table[a]+table[b]+1e-15) - (table[c] + 1e-15)/(table[c]+table[d]+1e-15))
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
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  double n = (double) M1.nrow();
  
  return (table[a]*table[d] - table[b]*table[c])/
    sqrt(n*(table[a] + table[b])*(table[a] + table[c]));
  
}

double syuleq(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return (table[a]*table[d] - table[b]*table[c])/(table[a]*table[d] + table[b]*table[c]);
  
}

double syuleqw(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return (sqrt(table[a]*table[d]) - sqrt(table[b]*table[c]))/
    (sqrt(table[a]*table[d]) + sqrt(table[b]*table[c]));
  
}

double dyuleq(
  const IntegerMatrix & M1,
  const IntegerMatrix & M2,
  bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return 2.0*table[b]*table[c]/(table[a]*table[d] + table[b]*table[c]);
  
}

//' @name similarity
//' @rdname similarity
//' @details (68) S Michael
//' @aliases Michael
double smichael(
  const IntegerMatrix & M1,
  const IntegerMatrix & M2,
  bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return 4.0*(table[a]*table[d] - table[b]*table[c])/
    (pow(table[a] + table[d], 2.0) + pow(table[b] + table[c], 2.0));
}

//' @name similarity
//' @rdname similarity
//' @details (73) Peirce
//' @aliases Peirce
double speirce(
  const IntegerMatrix & M1,
  const IntegerMatrix & M2,
  bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return (table[a]*table[b] + table[b]*table[c])/
    (table[a]*table[b] + 2*table[b]*table[c] + table[c]*table[d]);
  
}

// -----------------------------------------------------------------------------

typedef double (*funcPtr)(const IntegerMatrix & M1, const IntegerMatrix & M2, bool normalize);

NumericMatrix allsimilarities(
    const ListOf<IntegerMatrix> & M,
    bool normalized = false,
    funcPtr fun = NULL
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
        ans.at(pos++) = fun(M[i], M[j], normalized);
      }
      
      return cbind(I,J,ans);
      
}

void getmetric(std::string statistic, funcPtr & fun) {
  
  if      (statistic == "S14")      fun = &S14;
  else if (statistic == "hamming")  fun = &hamming;
  else if (statistic == "dennis")   fun = &dennis;
  else if (statistic == "starwid")  fun = &starwid;
  else if (statistic == "syuleq")   fun = &syuleq;
  else if (statistic == "syuleqw")  fun = &syuleqw;
  else if (statistic == "dyuleq")   fun = &dyuleq;
  else if (statistic == "smichael") fun = &smichael;
  else if (statistic == "speirce")  fun = &speirce;
  else Rcpp::stop("The statistic '%s' is not defined.", statistic);
  
  return ;
}

// [[Rcpp::export(name=".similarity")]]
NumericMatrix similarity(
    const ListOf<IntegerMatrix> & M,
    const std::string & statistic,
    bool normalized=false
  ) {
  
  funcPtr fun;
  getmetric(statistic, fun);
  
  return allsimilarities(M, normalized, fun);
  
}

// No longer needed
#undef a
#undef b
#undef c
#undef d
