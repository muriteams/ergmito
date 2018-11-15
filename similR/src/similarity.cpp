#include <Rcpp.h>
using namespace Rcpp;

#define a 0
#define b 1
#define c 2
#define d 3

typedef std::vector< int > vecint;

template<typename Ti, typename Tm> inline
void contingency_matrix(std::vector<Ti> & table, const Tm & M1, const Tm & M2) {
  
  // Catching error
  if ((M1.ncol() != M2.ncol()) | (M1.nrow() != M2.nrow())) 
    stop("`a` and `b` must have the same dimensions.");
  
  if (table.size() != 4)
    stop("`table` must be of size 4.");
  
  std::fill(table.begin(), table.end(), 0);
  
  int n = M1.nrow();
  int m = M1.ncol();
  
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      if (i == j) continue;
      else {
        if ((M1(i,j) == M2(i,j)) & (M1(i,j) == 1)) table[a]++;
        else if ((M1(i,j) != M2(i,j)) & (M1(i,j) == 1)) table[b]++;
        else if ((M1(i,j) != M2(i,j)) & (M1(i,j) == 0)) table[c]++;
        else if ((M1(i,j) == M2(i,j)) & (M1(i,j) == 0)) table[d]++;
      }
  
}


template<typename Ti, typename Tm> inline
std::vector<Ti> contingency_matrix(const Tm & M1, const Tm & M2) {
  
  std::vector<Ti> table(4);
  contingency_matrix< Ti, Tm >(table, M1, M2);
  
  return table;
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Jaccard (1): `"sjaccard"` or `"jaccard"`
//' @aliases Jaccard
double sjaccard(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalize = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return table[a]/(table[a] + table[b] + table[c]);
  
}


//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Faith (10): `"sfaith"` or `"faith"`
//' @aliases Faith
double sfaith(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalize = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return (table[a] + table[d]*0.5)/(table[a] + table[b] + table[c] + table[d]);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Gower & Legendre (11): `"sgl"` or `"gl"`
//' @aliases Gower-&-Legendre
double sgl(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalize = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return (table[a] + table[d])/(table[a] + 0.5*(table[b] + table[c]) + table[d]);
  
}

//' @name similarity
//' @rdname similarity
//' @section Distance:
//' - Sized Difference (24): `"dsd"` or `"sd"`
//' @aliases Sized-Difference
double dsd(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalize = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return pow(table[b] + table[c], 2.0)/
    pow(table[a] + table[b] + table[c] + table[d], 2.0);
  
}

//' @name similarity
//' @rdname similarity
//' @section Distance:
//' - Shaped Difference (25): `"dsphd"` or `"sphd"`
//' @aliases Shape-Difference
double dsphd(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalize = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return (((double) M1.nrow())*(table[b] + table[c]) - pow(table[b] - table[c], 2.0))/
    pow(table[a] + table[b] + table[c] + table[d], 2.0);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Tarwid (54): `"starwid"` or `"tarwid"`.
//' @aliases Tarwid
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

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Pearson & Heron 1 (54): `"sph1"` or `"ph1"` or `"s14"`. This is also known as S14 in
//'    Gower and Legendre (1986).
//'    
//'    In the case of the `S14` function, following Krackhardt's 1989:
//'    
//'    \deqn{%
//'    \sqrt{\left(\frac{a}{(a + c)} - \frac{b}{(b + d)}\right)\times\left(\frac{a}{(a + b)} - \frac{c}{(c + d)}\right)}
//'    }{%
//'    S14 = [(a/(a + c) - b/(b + d))*(a/(a + b) - c/(c + d))]^(1/2)
//'    }
//'   
//'    Which is an statistic lying between 0 and 1.
//'  
//' @aliases Person-&-Heron
//' @aliases S14
double sph1(
    const IntegerMatrix & R1,
    const IntegerMatrix & R2,
    bool normalized = false
  ) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(R1, R2);
  
  return (table[a]*table[d] - table[b]*table[c])/
    sqrt((table[a] + table[b])*(table[a] + table[c])*(table[b] + table[d])*(table[c] + table[d]));
        
}

//' @name similarity
//' @rdname similarity
//' @section Distance:
//' - Hamming (15): `"dhamming"` or `"hamming"`
//' @aliases Hamming
double dhamming(
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

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Mean Manhattan (20): `"dmh"` or `"mh"`
//'    \deqn{%
//'      D_{Mean-manhattan} = \frac{b + c}{a + b + c + d}
//'    }{%
//'     dmh = (b + c)/(a + b + c + d)
//'    }
//' @aliases Mean-Manhattan
double dmh(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return (table[b]+table[c])/
    (table[a] + table[b] + table[c] + table[d]);
    
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Dennis (44): `"sdennis"` or `"dennis"`
//' @aliases Dennis
double sdennis(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  double n = (double) M1.nrow();
  
  return (table[a]*table[d] - table[b]*table[c])/
    sqrt(n*(table[a] + table[b])*(table[a] + table[c]));
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Yuleq (61): `"syuleq"`
//' @aliases Yuleq
double syuleq(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return (table[a]*table[d] - table[b]*table[c])/(table[a]*table[d] + table[b]*table[c]);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Yuleq (63): `"syuleqw"`
//' @aliases Yuleq
double syuleqw(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return (sqrt(table[a]*table[d]) - sqrt(table[b]*table[c]))/
    (sqrt(table[a]*table[d]) + sqrt(table[b]*table[c]));
  
}

//' @name similarity
//' @rdname similarity
//' @section Distance:
//' - Yuleq (62): `"dyuleq"`
//' @aliases Yuleq
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
//' @section Similarity:
//' -  Michael (68):  `"smichael"` or `"michael"`
//'    
//'    \deqn{%
//'      S_{michael} = \frac{4(ad-bc)}{(a+d)^2 + (b+c)^2}
//'    }{%
//'      smichael = 4*(a*d - b*c)/[(a + d)^2 + (b + c)^2]
//'    }
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
//' @section Similarity:
//' -  Dispersion (66):  `"sdisp"` or `"disp"`
//'    
//'    \deqn{%
//'      S_{Dispersion} = \frac{ad - bc}{(a + b + c + d)^2}
//'    }{%
//'      sdisp = [a*d - b*c]/(a + b + c + d)^2
//'    }
//' @aliases Dispersion
double sdisp(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return (table[a]*table[d] - table[b] * table[c])/
    pow(table[a] + table[b] + table[c] + table[d], 2.0);
  
}


//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Hamann (67):  `"shamann"` or `"hamann"`
//'    
//'    \deqn{%
//'      S_{Hamann} = \frac{(a + d) - (b + c)}{a + b + c + d}
//'    }{%
//'      shamann = [(a + d) - (b + c)](a + b + c + d)
//'    }
//' @aliases Hamann
double shamann(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  
  return ((table[a] + table[d]) - (table[b] - table[c]))/
    (table[a] + table[b] + table[c] + table[d]);
  
}

// Auxiliar functions for Goodman & Kruskal, and Anderberg
inline double sigma(const std::vector<double> & table) {
  
  return std::max(table[a], table[b]) + std::max(table[c], table[d]) + 
    std::max(table[a], table[c]) + std::max(table[b], table[d]);
  
}

inline double sigma_prime(const std::vector<double> & table) {
  
  return std::max(table[a] + table[c], table[b] + table[d]) +
    std::max(table[a] + table[b], table[c] + table[d]);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Goodman & Kruskal (69): `"sgk"` or `"gk"`
//'    
//'    \deqn{%
//'      S_{Goodman&Kruskal} = \frac{\sigma - \sigma'}{2n - \sigma'} 
//'    }{%
//'      sgk = (s + s')/(2*n - s') 
//'    }
//'    
//'    where
//'    \eqn{\sigma = \max(a,b) + \max(c,d) + \max(a,c) + \max(b,d)}{
//'    s = max(a,b) + max(c,d) + max(a,c) + max(b,d)
//'    }, and 
//'    \eqn{\sigma' = \max(a + c, b + d) + \max(a + b, c + d)}{
//'    s' = max(a + c, b + d) + max(a + b, c + d)
//'    }
//'   
//' @aliases Goodman-&-Kruskal
double sgk(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  double s = sigma(table);
  double s_prime = sigma_prime(table);
  
  return (s - s_prime)/(2.0*(double) M1.nrow() - s_prime);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Anderberg (70): `"sanderberg"` or `"anderberg"`
//'    
//'    \deqn{%
//'      S_{Anderberg} = \frac{\sigma - \sigma'}{2n} 
//'    }{%
//'      sanderberg = (s + s')/(2*n) 
//'    }
//'    
//'    where \eqn{\sigma}{s} and \eqn{\sigma}{s'} are defined as in (69).
//'   
//' @aliases Anderberg
double sanderberg(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool normalized = false
) {
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(M1, M2);
  double s = sigma(table);
  double s_prime = sigma_prime(table);
  
  return (s - s_prime)/(2.0*(double) M1.nrow());
  
}





//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Peirce (73): `"speirce"` or `"peirce"`
//'    
//'    \deqn{%
//'      S_{Peirce} = \frac{ab + bc}{ab + 2bc + cd}
//'    }{%
//'      speirce = (a*b + b*c)/(a*b + 2*b*c + c*d)
//'    }
//'   
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

// [[Rcpp::export]]
IntegerMatrix reduce_dim(IntegerMatrix & x, int k) {
  
  // No reduction, return the same
  if (k == -1)
    return x;
  
  // Preparing reduced version of the matrix
  IntegerMatrix ans(x.nrow() - 1, x.ncol() - 1);
  ans.fill(0);
  
  int i0=0, j0;
  for (int i = 0; i < x.nrow(); ++i) {
    
    if (i == k)
      continue;

    j0 = 0;    
    for (int j = 0; j < x.ncol(); ++j) {
      
      // Self
      if (j == k) 
        continue;
      
      if (x(i,j) == 1)
        ans(i0, j0) = 1;

      // Out of range (can't do more)      
      if (++j0 == ans.ncol())
        break;
     
    }
    
    // Out of range (can't do more)
    if (++i0 == ans.nrow())
      break;
  }
  
  return ans;
  
}

// Applies whatever similarity/distance metric should be applying to all the
// requested combinations
NumericMatrix allsimilarities(
    const ListOf<IntegerMatrix> & M,
    bool normalized = false,
    funcPtr fun     = NULL,
    bool firstonly  = false,
    bool exclude_j    = false
) {
  
  int N  = M.size();
  int NN = firstonly? N-1 : N*(N-1)/2;
  NumericVector ans(NN);
  NumericVector I(NN),J(NN);
  
  int pos = 0;

  int firstloop = firstonly ? 1 : N;  
  if (exclude_j) {
    
    for (int i = 0; i < firstloop; ++i) 
      for (int j = i; j < N; ++j) {
        
        // Getting the pointers
        IntegerMatrix Mi = M[i];
        IntegerMatrix Mj = M[j];
        
        if (i == j) continue;
        else {
          I[pos] = i + 1;
          J[pos] = j + 1;
          ans[pos++] =
            fun(reduce_dim(Mi, j), reduce_dim(Mj, j), normalized);
        }
        
      }
        
  } else {
    
    for (int i = 0; i < firstloop; ++i) {
      for (int j = i; j < N; ++j)
        if (i == j) continue;
        else {
          I[pos] = i + 1;
          J[pos] = j + 1;
          ans[pos++] = fun(M[i], M[j], normalized);
        }
    }
    
  }
      
  return cbind(I,J,ans);
      
}

void getmetric(std::string s, funcPtr & fun) {
  
  if      (s == "sph1" | s == "ph1" | s == "s14") fun = &sph1;
  else if (s == "dhamming" | s == "hamming")      fun = &dhamming;
  else if (s == "dennis" | s == "sdennis")        fun = &sdennis;
  else if (s == "starwid" | s == "tarwid")        fun = &starwid;
  else if (s == "syuleq")                         fun = &syuleq;
  else if (s == "syuleqw")                        fun = &syuleqw;
  else if (s == "dyuleq")                         fun = &dyuleq;
  else if (s == "smichael" | s == "michael")      fun = &smichael;
  else if (s == "speirce" | s == "peirce")        fun = &speirce;
  else if (s == "sjaccard" | s == "jaccard")      fun = &sjaccard;
  else if (s == "sgk" | s == "gyk")               fun = &sgk;
  else if (s == "sanderberg" | s == "anderberg")  fun = &sanderberg;
  else if (s == "shamann" | s == "hamann")        fun = &shamann;
  else if (s == "dmh" | s == "mh")                fun = &dmh;
  else if (s == "sfaith" | s == "faith")          fun = &sfaith;
  else if (s == "sgl" | s == "gl")                fun = &sgl;
  else if (s == "dsd" | s == "sd")                fun = &dsd;
  else if (s == "dsphd" | s == "sphd")            fun = &dsphd;
  else if (s == "sdisp" | s == "disp")            fun = &sdisp;
  else Rcpp::stop("The statistic '%s' is not defined.", s);
  
  return ;
}

// [[Rcpp::export(name=".similarity")]]
NumericMatrix similarity(
    const ListOf<IntegerMatrix> & M,
    const std::string & statistic,
    bool normalized = false,
    bool firstonly  = false,
    bool exclude_j  = false
  ) {
  
  funcPtr fun;
  getmetric(statistic, fun);
  
  return allsimilarities(M, normalized, fun, firstonly, exclude_j);
  
}

// No longer needed
#undef a
#undef b
#undef c
#undef d
