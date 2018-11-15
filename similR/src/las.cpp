#include <Rcpp.h>
using namespace Rcpp;

//' @name las
//' @rdname las
//' @section Rule:
//' When `rule = "interesection"` or `rule = "i"`, the intersection method is used.
//' A tie is marked as present iff both i and j agree on the existance if (i,j).
void las_intersection(
    IntegerMatrix & LAS,
    const ListOf<IntegerMatrix> & M,
    int i,
    int j,
    int k = 0
) {
  
  // print(M);
  // print(LAS);
  // Rprintf("(i%i,j%i)\n",i,j);
  
  LAS(i, j) = (int) M[i](i, j)*M[j](i, j);
  
}

//' @name las
//' @rdname las
//' @section Rule:
//' When `rule = "union"` or `rule = "u"`, the union method is used. In such case
//' a tie is marked as present if any of the individuals declares it.
void las_union(
    IntegerMatrix & LAS,
    const ListOf<IntegerMatrix> & M,
    int i,
    int j,
    int k = 0
) {
  
  LAS(i, j) = (int) (M[i](i, j) | M[j](i, j));
  
}

//' @name las
//' @rdname las
//' @section Rule:
//' When `rule = "threshold"` or `rule = "t"`, the threshold method is used in which
//' a tie is marked as present iff the proportion of individuals who agree on the
//' existance of that is greater than `threshold`. In its paper, Krackhardt calls
//' this *Concensus Structure* (CS).
void las_threshold(
    IntegerMatrix & LAS,
    const ListOf<IntegerMatrix> & M,
    int i,
    int j,
    int k
) {
  
  LAS(i, j) += M[k](i, j);
  
}

typedef void (*lasFncPtr)(IntegerMatrix & LAS, const ListOf<IntegerMatrix> & M, int i, int j, int k);

void get_las(lasFncPtr & f, std::string rule) {
  
  if      (rule == "intersection" | rule == "i") f = &las_intersection;
  else if (rule == "union" | rule == "u") f = &las_union;
  else if (rule == "threshold" | rule == "t") f = &las_threshold;
  else
    stop("The LAS method '%s' is unkown.", rule);
  
}

// [[Rcpp::export(name=".las")]]
IntegerMatrix las(
    const ListOf< IntegerMatrix > & M,
    std::string rule,
    double threshold = 1,
    bool self = false
    ) {
  
  int nrow = M[0].nrow();
  int ncol = M[0].ncol();
  
  // Basic checks
  if (nrow != ncol)
    stop("LAS matrix can only be built in square matrices.");
  
  for (int i = 1; i < M.size(); ++i) {
    
    if (ncol != M[i].ncol())
      stop("Matrix %i has a different number of columns than matrix 1.",
           i + 1);
    
    if (nrow != M[i].nrow())
      stop("Matrix %i has a different number of rows than matrix 1.",
           i + 1);
    
  }
  
  if (M.size() != nrow)
    stop("There are not enough matrices to build the LAS.");
  
  // Getting the function
  lasFncPtr f;
  get_las(f, rule);
  
  IntegerMatrix LAS(nrow, ncol);
  LAS.fill(0);
  
  // Building the LAS
  if (rule != "t" && rule != "threshold") {
    
    for (int i = 0; i < nrow; ++i)
      for (int j = 0; j < nrow; ++j) {
        if (!self && i == j)
          continue;
        else
          f(LAS, M, i, j, 0);
      }
      
  } else {
    
    for (int j = 0; j < nrow; ++j)
      for (int i = 0; i < nrow; ++i)
        for (int k = 0; k < nrow; ++k) {
          if (!self && i == j)
            continue;
          else
            f(LAS, M, i, j, k);
        }
        
  }
      
  // If it was threshold based, then use it.
  if (rule == "threshold" | rule == "t")
    for (int i = 0; i < nrow; ++i)
      for (int j = 0; j < nrow; ++j) {
        
        if (!self && i == j) 
          continue;
        else 
          LAS(i,j) = (((double) LAS(i,j)/ (double) nrow) > threshold? 1 : 0);
        
      }
  
  return LAS;
  
}
