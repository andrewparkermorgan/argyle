#include <Rcpp.h>
using namespace Rcpp;

// compute proportion of alleles shared IBS between samples
// [[Rcpp::export]]
NumericMatrix dist_ibs(NumericMatrix x) {
   
   int nr = x.nrow();
   int nc = x.ncol();
   NumericMatrix D = NumericMatrix(nc, nc);
   
   for (int i = 0; i < nc; i++) { // indiv 1
    // check for user interrupt once every row
   	// of final distance matrix, not too expensive
   	checkUserInterrupt();
   	for (int j = i+1; j < nc; j++) { // indiv 2
   		NumericVector nna = ifelse( is_na(x(_,i)) | is_na(x(_,j)), 0.0, 1.0 );
   		D(j,i) = sum( abs(na_omit(x(_,j) - x(_,i))) )/(2*sum(nna));
   	}
   }
   
   return D;

}
