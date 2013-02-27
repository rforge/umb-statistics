#include <Rcpp.h>

RcppExport SEXP qnw_sg_mult(SEXP in1, SEXP in2, SEXP in3) {
    //BEGIN_RCPP

    Rcpp::List seqs1(in1);
    Rcpp::List seqs2(in2);
	Rcpp::IntegerVector my(in3);
    
    int ns1 = seqs1.length();
    int ns2 = seqs2.length();
    Rcpp::NumericMatrix dist(ns1,ns2);
    Rcpp::NumericMatrix dist_scaled(ns1,ns2);
    Rcpp::NumericMatrix J(ns1,ns2);
    std::vector<int> D;
	D.resize( my(0)+1 );
	std::vector<int> L;
	L.resize( my(0)+1 );
	int diagonal,horizontal,vertical,nx,ny,d,jj,ll,prevD,prevL,s;
	std::vector<int> seq1;
	std::vector<int> seq2;

    for( int k=0; k<ns1; k++) {
		seq1 = seqs1[k];
		nx = seq1.size();
		
        for( int l=0; l<ns2; l++) {
			seq2 = seqs2[l];
			ny = seq2.size();

			for( int j = 0; j <= ny; j++ ) D[ j ] = 0;
			for( int j = 0; j <= ny; j++ ) L[ j ] = 0;
			for( int i = 1; i <= nx; i++) {
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

			d = D[0];
			ll = L[0];
            jj = 1;
            for( int i = 1; i <= ny; i++ )
                if( d < D[i] ) {
                    d = D[i];
                    ll = L[i];
                    jj = i;
                }
            dist( k , l ) = d;
            J( k , l ) = jj;
			dist_scaled( k , l ) = 0.5 - d/(2.0*ll);
        }
    }
    return Rcpp::List::create(Rcpp::Named("D") = dist, Rcpp::Named("D_scaled") = dist_scaled, Rcpp::Named("j") = J); 

    //END_RCPP
}
