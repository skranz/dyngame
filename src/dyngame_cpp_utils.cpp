
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List C_get_from_pos(Environment m, IntegerMatrix pos_k) {
	int nx = m["nx"];
	int n = m["n"];
	
	//IntegerMatrix pos_k(nx,n+1);
	
	List opt = m["opt"];
	
	// Current pos in list for state x (row) state k (column)
	IntegerMatrix pos(nx,n+1);
	NumericVector G(nx),C(nx);
	NumericMatrix Lk(nx,n+1);
	IntegerMatrix ak(nx,n+1);
		
	for (int x=0;x<nx;x++) {
		List optx = opt[x];
		NumericMatrix mat=optx["mat.e"];
		int k = 0;
		Lk(x,k) = mat(pos(x,k),0);
		G(x) = mat(pos(x,k),1);
		ak(x,k) = mat(pos(x,k),2);

		double C_temp = 0;	
		for (int i=0; i<n;i++) {
			List mat_i = optx["mat.i"];
			NumericMatrix mat=mat_i[i];
			int k = i+1;
			Lk(x,k) = mat(pos(x,k),0);
			C_temp += mat(pos(x,k),1);
			ak(x,k) = mat(pos(x,k),2);
		}
		C[x] = C_temp;		
	}
	List ret = List::create(Named("G")=G,Rcpp::Named("C") = C,
		Named("ak")=ak,Rcpp::Named("Lk") = Lk);
	return ret;
}