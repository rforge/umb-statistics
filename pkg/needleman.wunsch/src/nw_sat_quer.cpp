#include <Rcpp.h>

RcppExport SEXP nw_sat_quer(SEXP in1, SEXP in2, SEXP in3, SEXP in4, SEXP in5, SEXP in6) {
    //BEGIN_RCPP

    Rcpp::List seqs1(in1);
    Rcpp::List seqs2(in2);
    Rcpp::IntegerVector s(in3);
    Rcpp::IntegerVector gap(in4);
	Rcpp::IntegerVector mx(in5);
	Rcpp::IntegerVector my(in6);
    
    int ns1 = seqs1.length();
    int ns2 = seqs2.length();
	int diagonal, horizontal, vertical;
	
	std::vector< std::vector<int> > D;
    D.resize( mx(0)+1 , std::vector<int>( my(0)+1 , 0 ) );
	std::vector< std::vector<int> > L;
    L.resize( mx(0)+1 , std::vector<int>( my(0)+1 , 0 ) );
	std::vector< std::vector<int> > S;
	S.resize( 4 , std::vector<int>( 4 , 0 ) );
    for( int i = 0; i<4; i++ ) {
        for( int j = 0; j<4; j++) {
            S[i][j] = s[i+j*4];
        }
    }
    Rcpp::NumericMatrix dist(ns1,ns2);
    Rcpp::NumericMatrix dist_scaled(ns1,ns2);

    for( int k=0; k<ns1; k++) {
        SEXP s1 = seqs1[k];
        Rcpp::IntegerVector seq1(s1);
        int nx = seq1.length();

        for( int l=0; l<ns2; l++) {
            SEXP s2 = seqs2[l];
            Rcpp::IntegerVector seq2(s2);
            int ny = seq2.length();
            
            for( int i = 0; i <= nx; i++ ) {
				for( int j = 0; j <= ny; j++ ) D[ i ][ j ] = 0;
			}
            for( int i = 0; i <= ny; i++ ) D[ 0 ][ i ] = i * gap[0];
            for( int i = 0; i <= nx; i++ ) D[ i ][ 0 ] = i * gap[0];
            for( int i = 0; i <= nx; i++ ) {
				for( int j = 0; j <= ny; j++ ) L[ i ][ j ] = 0;
			}
            for( int i = 0; i <= ny; i++ ) L[ 0 ][ i ] = i;
            for( int i = 0; i <= nx; i++ ) L[ i ][ 0 ] = i;
            
            for( int i = 1; i <= nx; i++) {
                for( int j = 1; j <= ny; j++ ){
                    diagonal = D[i-1][j-1] + S[seq1[i-1]][seq2[j-1]];
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
            dist( k , l ) = D[nx][ny];
			dist_scaled( k , l ) = 0.5 - D[nx][ny]/(2.0*L[nx][ny]);
//			minD = nx;
//			if( ny < nx ){
//				minD = ny;
//			}
//           dist_scaled( k , l ) = 0.5 * (1.0 - (D[nx][ny] + abs(nx-ny)) / minD);
        }
    }
    return Rcpp::List::create(Rcpp::Named("D") = dist, Rcpp::Named("D_scaled") = dist_scaled); 

    //END_RCPP
}
