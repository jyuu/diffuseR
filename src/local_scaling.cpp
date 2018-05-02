
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix scaledist(const NumericMatrix & x, NumericVector scale) {
	unsigned int outrows = x.nrow(), i = 0, j = 0; 
	double d; 
	NumericMatrix out(outrows, outrows); 

	for (i = 0; i < outrows - 1; i++){
		out(i,i) = 1;
		for (j = i+1; j < outrows; j++){
			d = x(i,j)/(scale(i)*scale(j));
			out(j, i) = exp(-d); 
			out(i, j) = exp(-d);
		}
	}
	
    return out;
}