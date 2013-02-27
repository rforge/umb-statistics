#include <Rcpp.h>

RcppExport SEXP qnw(SEXP in1, SEXP in2) {
    //BEGIN_RCPP

    Rcpp::IntegerVector seq1(in1);
    Rcpp::IntegerVector seq2(in2);

    int nx = seq1.length();
    int ny = seq2.length();
	int diagonal, horizontal, vertical, prevD, prevL, s, l;

    std::vector<int> D;
	D.resize( ny+1 );
    for( int j = 0; j <= ny; j++ ) D[ j ] = -1 * j;
	std::vector<int> L;
	L.resize( ny+1 );
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
			l = L[j-1] + 1;
            if( horizontal > diagonal ){
                diagonal = horizontal;
				l = prevL + 1;
			}
            if( vertical > diagonal ){
                diagonal = vertical;
				l = L[j] + 1;
			}
			D[j-1] = prevD;
			L[j-1] = prevL;
			prevD = diagonal;
			prevL = l;
		}
		D[ny] = prevD;
		L[ny] = prevL;
	}
    return Rcpp::List::create(Rcpp::Named("D") = D[ny], Rcpp::Named("L") = L[ny]);

    //END_RCPP
}
