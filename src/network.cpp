#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List init_network(
    int n,
    bool directed  = true,
    bool hyper     = false,
    bool loops     = false,
    bool multiple  = false,
    bool bipartite = false
) {
  
  if (n < 0)
    stop("`n` cannot be less than 0.");
  
  List emptylist;

  List g = List::create(
    /* Master edgelist:
     * This is in itself a list of length (mnext - 1) which, if different from
     * an empty list, has the following elements: inl (tail), outl (head), alt (attributes)
     * 
     * See https://github.com/statnet/network/blob/525286ff257745011f35a58c193c68745f78e853/vignettes/networkVignette.Rnw#L331
     */ 
    _["mel"] = List(0),
    
    _["gal"] = List::create(     // Graph attribute list
      _["n"]         = n,
      _["mnext"]     = 1,        // Next number of edges (next m)
      _["directed"]  = directed,
      _["hyper"]     = hyper,
      _["loops"]     = loops,
      _["multiple"]  = multiple,
      _["bipartite"] = bipartite
    ),
    _["val"]       = R_NilValue,
    _["iel"]       = R_NilValue,
    _["oel"]       = R_NilValue
    );
  
  // Depending if there are nodes, then
  if (n > 0) {

    std::vector< List > listoflists(n);
    std::vector< IntegerVector > listofInt(n);
    
    // n lists
    g["val"] = listoflists;
    g["iel"] = listofInt;
    g["oel"] = listofInt;


  } else {

    // n lists
    g["val"] = emptylist;
    g["iel"] = emptylist;
    g["oel"] = emptylist;

  }
  
  g.attr("class") = "network";
  
  return g;
  
}

inline ListOf< List > matrix_to_networki(
    const IntegerMatrix & x,
    bool directed  = true,
    bool hyper     = false,
    bool loops     = false,
    bool multiple  = false,
    bool bipartite = false
    ) {
  
  int n = x.nrow(), m = x.ncol();
  
  if (n != m)
    stop("Non-square matrix.");
  
  std::vector< List > ans;
  ans.reserve(x.size());
  
  std::vector< std::vector<int> > iel(n), oel(n);
  
  // Default lists. We do this so that we don't make that many calls during the
  // for-loop.
  List NA_trueList  = List::create(_["na"] = true);
  List NA_falseList = List::create(_["na"] = false);
    
  
  int nedge = 0;
  for (int j = 0; j < m; ++j) {
    
    int ni = n;
    if (!directed)
      ni = j;
    
    for (int i = 0; i < ni; ++i) {
      
      if (!loops && (i == j)) 
        continue;
        
      if (x(i, j) != 0) {
        
        List edge = List::create(
          _["inl"]  = i + 1, //
          _["outl"] = j + 1, // While we see from 0 to (n-1), R lives in 1 to n.
          _["atl"]  = ( x(i, j) == R_NaInt ) ? NA_trueList : NA_falseList
        );
        
        ans.push_back(edge);
        iel.at(j).push_back(++nedge);
        oel.at(i).push_back(nedge);
        
#ifdef ERGMITO_DEBUG_NETWORK
        Rprintf("x[%d, %d]: edge %d\n", i, j, nedge);
#endif
        
      }
    }
  }
      
  // typedef std::vector< std::vector< int > > vecvecint;
  // for (vecvecint::iterator iter = iel.begin(); iter != iel.end(); ++iter) {
  //   std::reverse(iter->begin(), iter->end());
  // }
  // for (vecvecint::iterator iter = oel.begin(); iter != oel.end(); ++iter) {
  //   std::reverse(iter->begin(), iter->end());
  // }
  
  // Empty list of attributes
  ListOf< List > listoflist(n);
  
  List network = List::create(
    _["mel"] = wrap(ans),
    _["gal"] = List::create(
      _["n"]         = n,
      _["mnext"]     = ++nedge,        // Next number of edges (next m)
      _["directed"]  = directed,
      _["hyper"]     = hyper,
      _["loops"]     = loops,
      _["multiple"]  = multiple,
      _["bipartite"] = bipartite
    ),
    _["val"] = listoflist,
    _["iel"] = wrap(iel),
    _["oel"] = wrap(oel)
  );
    
  network.attr("class") = "network";
  
  return network;
  
}

// [[Rcpp::export(name="matrix_to_network.")]]
ListOf< List > matrix_to_network(
    const ListOf< IntegerMatrix > & x,
    const LogicalVector & directed ,
    const LogicalVector & hyper ,
    const LogicalVector & loops ,
    const LogicalVector & multiple ,
    const LogicalVector & bipartite
  ) {
  
  int n = x.size();
  ListOf< List > out(n);
  
  for (int i = 0; i < n; ++i)
    out[i] = matrix_to_networki(
      x[i], directed[i], hyper[i], loops[i], multiple[i], bipartite[i]
  );
  
  return out;
  
}
