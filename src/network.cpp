#include <Rcpp.h>
using namespace Rcpp;

// #define ERGMITO_DEBUG_NETWORK 1

// // [[Rcpp::export]]
// List init_network(
//     int n,
//     bool directed  = true,
//     bool hyper     = false,
//     bool loops     = false,
//     bool multiple  = false,
//     bool bipartite = false
// ) {
//   
//   if (n < 0)
//     stop("`n` cannot be less than 0.");
//   
//   List emptylist;
// 
//   List g = List::create(
//     /* Master edgelist:
//      * This is in itself a list of length (mnext - 1) which, if different from
//      * an empty list, has the following elements: inl (tail), outl (head), alt (attributes)
//      * 
//      * See https://github.com/statnet/network/blob/525286ff257745011f35a58c193c68745f78e853/vignettes/networkVignette.Rnw#L331
//      */ 
//     _["mel"] = List(0),
//     
//     _["gal"] = List::create(     // Graph attribute list
//       _["n"]         = n,
//       _["mnext"]     = 1,        // Next number of edges (next m)
//       _["directed"]  = directed,
//       _["hyper"]     = hyper,
//       _["loops"]     = loops,
//       _["multiple"]  = multiple,
//       _["bipartite"] = bipartite
//     ),
//     _["val"]       = R_NilValue,
//     _["iel"]       = R_NilValue,
//     _["oel"]       = R_NilValue
//     );
//   
//   // Depending if there are nodes, then
//   if (n > 0) {
// 
//     std::vector< List > listoflists(n);
//     std::vector< IntegerVector > listofInt(n);
//     
//     // n lists
//     g["val"] = listoflists;
//     g["iel"] = listofInt;
//     g["oel"] = listofInt;
// 
// 
//   } else {
// 
//     // n lists
//     g["val"] = emptylist;
//     g["iel"] = emptylist;
//     g["oel"] = emptylist;
// 
//   }
//   
//   g.attr("class") = "network";
//   
//   return g;
//   
// }

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
  
  std::vector< List > mel, val(n);
  mel.reserve(x.size());
  
  std::vector< std::vector<int> > iel(n), oel(n);
  
  // Default lists. We do this so that we don't make that many calls during the
  // for-loop.
  List NA_trueList  = List::create(_["na"] = true);
  List NA_falseList = List::create(_["na"] = false);
  
  int nedge = 0;
  for (int j = 0; j < m; ++j) {
    
    char name[64];
    sprintf(&(name[0]), "%i", (unsigned short) j + 1u);
    
    val.at(j) = NA_falseList;
    val.at(j)["vertex.names"] = name;
    
    int ni = n;
    if (!directed)
      ni = j;
    
    for (int i = 0; i < ni; ++i) {
      
      if (!loops && (i == j)) 
        continue;
        
      if (x(i, j) != 0) {
        
        List edge = List::create(
          _["inl"]  = j + 1, //
          _["outl"] = i + 1, // While we see from 0 to (n-1), R lives in 1 to n.
          _["atl"]  = ( x(i, j) == R_NaInt ) ? NA_trueList : NA_falseList
        );
        
        mel.push_back(edge);
        iel.at(j).push_back(++nedge);
        oel.at(i).push_back(nedge);
        
#ifdef ERGMITO_DEBUG_NETWORK
        Rprintf("[matrix_to_network] x[%d, %d]: edge %d\n", i, j, nedge);
#endif
        
      }
    }
  }
  
  // Resizing back
  mel.shrink_to_fit();
      
  // typedef std::vector< std::vector< int > > vecvecint;
  // for (vecvecint::iterator iter = iel.begin(); iter != iel.end(); ++iter) {
  //   std::reverse(iter->begin(), iter->end());
  // }
  // for (vecvecint::iterator iter = oel.begin(); iter != oel.end(); ++iter) {
  //   std::reverse(iter->begin(), iter->end());
  // }
  
  List network = List::create(
    _["mel"] = mel,
    _["gal"] = List::create(
      _["n"]         = n,
      _["mnext"]     = ++nedge,        // Next number of edges (next m)
      _["directed"]  = directed,
      _["hyper"]     = hyper,
      _["loops"]     = loops,
      _["multiple"]  = multiple,
      _["bipartite"] = bipartite
    ),
    _["val"] = val,
    _["iel"] = iel,
    _["oel"] = oel
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
// 
// inline ListOf< List > add_vertex_attri(
//     const ListOf< List > & x, 
//     const ListOf< GenericVector > & vattrs,
//     const StringVector & names,
//     const List & newattrs
//   ) {
//   
//   ListOf< List > val = x;
//   
// #ifdef ERGMITO_DEBUG_NETWORK
//   Rprintf("[add_vertex_attr] has class?: %d\n", (int) as<List>(x).hasAttribute("class"));
// #endif
//   
//   // Checking lengths
//   int nattrs = vattrs.size(), n = val.size();
//   for (int i = 0; i < nattrs; ++i)
//     if (vattrs[i].size() != n)
//       stop("The length of `vattr` doesn't matches the number of vertices in the graph.");
//   
//   int nattrs_old = val[0].size();
//   for (int i = 0; i < n; ++i) {
//     
//     // Appending the old attributes
//     List newattrs_i = clone(newattrs);
//     
//     for (int j = 0; j < nattrs_old; ++j)
//       newattrs_i[j] = val[i].at(j);
//     
//     // Adding the attributes to the ith vertex
//     for (int j = 0; j < nattrs; ++j) {
//       newattrs_i[j + nattrs_old] = vattrs[j].at(i);
//     }
//     
// #ifdef ERGMITO_DEBUG_NETWORK
//     Rprintf("[add_vertex_attr] This is how the new attributes look like of i=%d:\n", i);
//     print(wrap(newattrs_i));
// #endif
//     
//     // We need to update the names
//     val[i] = newattrs_i;
//     
//   }
//   
//   // Updating the vertex attributes in x
//   return val;
//   
// }
// 
// // [[Rcpp::export(name="add_vertex_attr.")]]
// ListOf< List > add_vertex_attr(
//     const ListOf< List > & x, 
//     const ListOf< GenericVector > & vattrs,
//     const StringVector & names,
//     const ListOf< List > & oldattrs
//   
// ) {
//   
//   int n = x.size();
//   ListOf< List > res(n);
//   
//   // Checking the length
//   if (vattrs.size() != names.size())
//     stop("The number of names must match the length of a.");
//   
//     
// #ifdef ERGMITO_DEBUG_NETWORK
//     Rprintf("[add_vertex_attr] Going OK up to basic checks. Found %d attributes to add.\n", names.size());
// #endif
//     
//     // Capturing names
//     StringVector newnames = oldattrs.names();
//     
// #ifdef ERGMITO_DEBUG_NETWORK
//     Rprintf("[add_vertex_attr] Here is the lit of attributes to be added:\n");
//     print(wrap(newnames));
// #endif
//     
//     for (StringVector::const_iterator iter = names.begin(); iter != names.end(); ++iter)
//       newnames.push_back(*iter);
// 
// #ifdef ERGMITO_DEBUG_NETWORK
//     Rprintf("[add_vertex_attr] Creating the `newnames` object went fine.\n");
// #endif
//         
//     // All will have the same number of attributes
//     List newattrs(newnames.size());
//     newattrs.names() = newnames;
//     
// #ifdef ERGMITO_DEBUG_NETWORK
//     Rprintf("[add_vertex_attr] Starting with the loop.\n");
// #endif
//   
//   for (int i = 0; i < n; ++i) {
//     List x_i = clone(as< List >(x[i]));
//     // ListOf< List > val = as< ListOf< List > >(x_i["val"]);
//     // 
//     // // Updating the value
//     // add_vertex_attri(val, vattrs, names);
//     // x_i["val"] = val;
//     x_i["val"] = add_vertex_attri(
//       as< ListOf< List > >(x_i["val"]), vattrs, names, newattrs
//       );
// 
//     res[i] = x_i;
//   }
//   
//   return res;
//   
// }