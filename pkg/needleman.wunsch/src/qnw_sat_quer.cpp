#include <Rcpp.h>

RcppExport SEXP qnw_sat_quer(SEXP in1, SEXP in2, SEXP in3) {
    //BEGIN_RCPP

    Rcpp::List seqs1(in1);
    Rcpp::List seqs2(in2);
	Rcpp::IntegerVector my(in3);
    
    int ns1 = seqs1.length();
    int ns2 = seqs2.length();
    Rcpp::NumericMatrix dist(ns1,ns2);
    Rcpp::NumericMatrix dist_scaled(ns1,ns2);
    std::vector<int> D;
	D.resize( my(0)+1 );
	std::vector<int> L;
	L.resize( my(0)+1 );
	int s,diagonal,vertical,horizontal,nx,ny,prevD,prevL,ll;
	std::vector<int> seq1;
	std::vector<int> seq2;

	for( int k=0; k<ns1; k++) {
		seq1 = seqs1[k];
		nx = seq1.size();
		
		for( int l=0; l<ns2; l++) {
			seq2 = seqs2[l];
			ny = seq2.size();
			
			for( int j = 0; j <= ny; j++ ) D[ j ] = -1 * j;
			for( int j = 0; j <= ny; j++ ) L[ j ] = j;
			for( int i = 1; i <= nx; i++ ){
				prevD = -1 * i;
				prevL = i;
				for( int j = 1; j <= ny; j++ ){
					s = -1;
					if( seq1[i-1] == seq2[j-1] )
						s = 1;
					diagonal = D[j-1] + s;
					horizontal = prevD - 1;
					vertical = D[j] - 1;
					ll = L[j-1] + 1;
					if( horizontal > diagonal ){
						diagonal = horizontal;
						ll = prevL + 1;
					}
					if( vertical > diagonal ){
						diagonal = vertical;
						ll = L[j] + 1;
					}
					D[j-1] = prevD;
					L[j-1] = prevL;
					prevD = diagonal;
					prevL = ll;
				}
				D[ny] = prevD;
				L[ny] = prevL;
			}

            dist( k , l ) = D[ny];
			dist_scaled( k , l ) = 0.5 - D[ny]/(2.0*L[ny]);
        }
    }
    return Rcpp::List::create(Rcpp::Named("D") = dist, Rcpp::Named("D_scaled") = dist_scaled); 

    //END_RCPP
}
