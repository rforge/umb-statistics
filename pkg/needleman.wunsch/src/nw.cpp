#include <Rcpp.h>

RcppExport SEXP nw(SEXP in1, SEXP in2, SEXP in3, SEXP in4) {
    //BEGIN_RCPP

    Rcpp::IntegerVector seq1(in1);
    Rcpp::IntegerVector seq2(in2);
    Rcpp::IntegerVector s(in3);
    Rcpp::IntegerVector gap(in4);

    int nx = seq1.length();
    int ny = seq2.length();
	int diagonal, horizontal, vertical;
            

	std::vector< std::vector<int> > D;
	D.resize( nx+1 , std::vector<int>( ny+1 , 0 ) );
    for( int i = 0; i <= ny; i++ ) D[ 0 ][ i ] = i * gap[0];
    for( int i = 0; i <= nx; i++ ) D[ i ][ 0 ] = i * gap[0];
	std::vector< std::vector<int> > L;
	L.resize( nx+1 , std::vector<int>( ny+1 , 0 ) );
    for( int i = 0; i <= ny; i++ ) L[ 0 ][ i ] = i;
    for( int i = 0; i <= nx; i++ ) L[ i ][ 0 ] = i;

	std::vector< std::vector<int> > S;
	S.resize( 4 , std::vector<int>( 4 , 0 ) );
    for( int i = 0; i<4; i++ ) {
        for( int j = 0; j<4; j++) {
            S[i][j] = s[i+j*4];
        }
    }
    
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
    return Rcpp::List::create(Rcpp::Named("D") = D[nx][ny], Rcpp::Named("L") = L[nx][ny]);

    //END_RCPP
}
