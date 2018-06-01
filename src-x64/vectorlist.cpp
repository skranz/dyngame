#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector C_RowSums_VectorList(NumericVector vec, IntegerVector ncols) {
		NumericVector xvec(vec);                
		IntegerVector xncols(ncols);
		
		int nrows = xncols.length();
		
		NumericVector val(nrows);

		double* pvec = xvec.begin();
		int*    pncols = xncols.begin();
		double* pval = val.begin();

		for (int i = 0; i < nrows; i++) {
			double sum = 0;
			for (int j = 0; j < pncols[i]; j++) {
				sum += *pvec;
				++pvec;
			}
			pval[i] = sum;
		}
		return val;
}

// [[Rcpp::export]]
NumericVector C_RowMaxs_VectorList(NumericVector vec,IntegerVector ncols) {
		int nrows = ncols.length();	
		NumericVector val(nrows);

		double* pvec = vec.begin();
		int*    pncols = ncols.begin();
		double* pval = val.begin();

		for (int i = 0; i < nrows; i++) {
			if (pncols[i]>0) {
				double max = *pvec;
				++pvec;
				for (int j = 1; j < pncols[i]; j++) {
					if (*pvec>max) {
						max = *pvec;
					}
					++pvec;
				}
				pval[i] = max;
			} else {
				pval[i] = 0;
			}
		}
		return val;
}

// [[Rcpp::export]]
IntegerVector C_which_RowMaxs_VectorList(NumericVector vec,IntegerVector ncols) {
		int nrows = ncols.length();
		
		IntegerVector val(nrows);

		double* pvec = vec.begin();
		int*    pncols = ncols.begin();
		int*    pval = val.begin();
	 
		int counter = 0;
		for (int i = 0; i < nrows; i++) {
			
			if (pncols[i]>0) {
				counter++;
				double max = *pvec;
				int maxind = counter;
				pvec++;
				
				for (int j = 1; j < pncols[i]; j++) {
					counter++;
					if (*pvec > max) {
						max = *pvec;
						maxind = counter;
					} 
					pvec++;
				}
				pval[i] = maxind;
			} else {
				pval[i] = 0;
			}
			
		}
		
		return val;
}

	

