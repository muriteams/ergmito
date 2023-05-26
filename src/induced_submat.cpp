#include <Rcpp.h>

using namespace Rcpp;

inline IntegerMatrix induced_submati(
    const IntegerMatrix & net,
    const IntegerVector & v
) {
  
  if (net.nrow() != net.ncol())
    stop("`net` should be a square matrix");
  
  IntegerVector v_unique = Rcpp::unique(v);
  
  if (v.size() != v_unique.size())
    stop("`v` has repeated elements.");
  
  unsigned int i, j, n = v.size();
  IntegerMatrix newnet(n, n);
  newnet.fill(0);
  
  for (i = 0u; i < n; ++i) 
    for (j = 0u; j < n; ++j) {
      
      if ((v[i] < 0) || (v[i] > (net.size() - 1)))
        stop("Vertex index out of range");
      if ((v[j] < 0) || (v[j] > (net.size() - 1)))
        stop("Vertex index out of range");
        
      if (net(v[i], v[j]) != 0)
        newnet(i, j) = net(v[i], v[j]);
    }
      
  return newnet;
  
}

// [[Rcpp::export(name = "induced_submat_cpp")]]
std::vector< IntegerMatrix > induced_submat(
  const std::vector< IntegerMatrix > & nets,
  const std::vector< IntegerVector > & vs
) {
  
  // Initializing variables
  unsigned int * net_ptr, * vs_ptr, *N,
  i, zero = 0u, nnets = nets.size(), nvs = vs.size();
  
  if ((nnets == 0) | (nvs == 0)) {
    stop("One of `nets` or `vs` is of length 0.");
  } else  if (((nnets > 1) & (nvs > 1)) && (nets.size() != vs.size())) {
    stop("The length of `nets` should be the same as `vs`.");
  } else if (nnets == nvs) {
    N       = &nnets;
    net_ptr = &i;
    vs_ptr  = &i;
  } else if ((nnets > 1) & (nvs == 1)) {
    N       = &nnets;
    net_ptr = &i;
    vs_ptr  = &zero;
  } else if ((nnets == 1) & (nvs > 1)) {
    N       = &nvs;
    net_ptr = &zero;
    vs_ptr  = &i;
  } else
    stop("Both nets and vertices have more than one element and different sizes.");
  
  std::vector< IntegerMatrix > ans;
  ans.reserve(*N);
  
  for (i = 0u; i < *N; ++i)
    ans.push_back(induced_submati(nets[*net_ptr], vs[*vs_ptr] - 1));
  
  return ans;
  
}


  
