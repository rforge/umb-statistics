#include <Rcpp.h>

RcppExport SEXP nws_sg_mult(SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP in6) {
    //BEGIN_RCPP

    Rcpp::List seqs1(in1);
    Rcpp::List seqs2(in2);
    Rcpp::IntegerVector S(in3);
    Rcpp::IntegerVector gap(in4);
	Rcpp::IntegerVector mx(in5);
	Rcpp::IntegerVector my(in6);
    
    int ns1 = seqs1.length();
    int ns2 = seqs2.length();
	int diagonal, horizontal, vertical, nx, ny, d, jj, ll;
//	double tmp;
    Rcpp::NumericMatrix dist(ns1,ns2);
    Rcpp::NumericMatrix dist_scaled(ns1,ns2);
    Rcpp::NumericMatrix J(ns1,ns2);
	std::vector< std::vector<int> > D;
    D.resize( mx(0)+1 , std::vector<int>( my(0)+1 , 0 ) );
	std::vector< std::vector<int> > L;
    L.resize( mx(0)+1 , std::vector<int>( my(0)+1 , 0 ) );

    for( int k=0; k<ns1; k++) {
        for( int l=0; l<ns2; l++) {
            SEXP s1 = seqs1[k];
            Rcpp::IntegerVector seq1(s1);
            SEXP s2 = seqs2[l];
            Rcpp::IntegerVector seq2(s2);
            nx = seq1.length();
            ny = seq2.length();
            
            for( int i = 0; i <= nx; i++ ) {
				for( int j = 0; j <= ny; j++ ) D[ i ][ j ] = 0;
			}
            for( int i = 0; i <= nx; i++ ) D[ i ][ 0 ] = i * gap[0];
            for( int i = 0; i <= nx; i++ ) {
				for( int j = 0; j <= ny; j++ ) L[ i ][ j ] = 0;
			}
            for( int i = 0; i <= nx; i++ ) L[ i ][ 0 ] = i;
            
            for( int i = 1; i <= nx; i++) {
                for( int j = 1; j <= ny; j++ ){
                    int s = S[1];
                    if( seq1[i-1] == seq2[j-1] )
                        s = S[0];
                    diagonal = D[i-1][j-1] + s;
                    horizontal = D[i][j-1] + gap[0];
                    vertical = D[i-1][j] + gap[0];
					L[i][j] = L[i-1][j-1] + 1;
                    if( horizontal > diagonal ){
                        diagonal = horizontal;
						L[i][j] = L[i][j-1] + 1;
					}
                    if( vertical > diagonal ){
                        diagonal = vertical;
						L[i][j] = L[i-1][j] + 1;
					}
                    D[i][j] = diagonal;
                }
            }
            d = D[nx][0];
			ll = L[nx][0];
            jj = 1;
            for( int i = 1; i <= ny; i++ )
                if( d < D[nx][i] ) {
                    d = D[nx][i];
                    ll = L[nx][i];
                    jj = i;
                }
            dist( k , l ) = d;
            J( k , l ) = jj;
			dist_scaled( k , l ) = 0.5 - d/(2.0*ll);
//            tmp = (d + nx);
//            tmp = tmp / 2.0;
//            tmp = tmp / nx;
//            dist_scaled( k , l ) = 1.0 - tmp;
        }
    }
    return Rcpp::List::create(Rcpp::Named("D") = dist, Rcpp::Named("D_scaled") = dist_scaled, Rcpp::Named("j") = J); 

    //END_RCPP
}
