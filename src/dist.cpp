#include <Rcpp.h>
using namespace Rcpp;

// compute proportion of alleles shared IBS between samples
// [[Rcpp::export]]
NumericMatrix dist_ibs(NumericMatrix x) {
   
   int nr = x.nrow();
   int nc = x.ncol();
   NumericMatrix D = NumericMatrix(nc, nc);
   
   for (int i = 0; i < nc; i++) {
   	for (int j = i+1; j < nc; j++) {
   		D(j,i) = sum(abs(na_omit(x(_,i)-x(_,j)))) / (2*nr);
   	}
   }
   
   return D;

}
