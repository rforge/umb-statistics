#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector Kmer_classify_alien_RDP( SEXP Qin, SEXP seqs, int K ) {
  
  NumericMatrix Q(Qin);
  Rcpp::List strings(seqs);
  int nclass = Q.nrow();
  int where = 0;
  int nElem = pow(4,K);
  int num_strings = strings.length();
  NumericVector p(nElem);
  NumericVector Ci(nclass);
  LogicalVector X(nElem);
  IntegerVector C(num_strings);
  std::vector<int> Where(K);

  // Apply log2 to training set
  for(int i=0; i < Q.nrow(); i++){
    for(int j=0; j < Q.ncol(); j++){
      Q(i,j) = log2(Q(i,j));
    }
  }

  // Prepare powers of 4
  for(int i=0; i < K; i++){
    Where[i] = pow(4,K-i-1);
  }
  
  // Loop over sequences
  for( int i=0; i < num_strings; i++ ) {
    SEXP s1 = strings[i];
    Rcpp::IntegerVector seq1(s1);
    int num_substr = seq1.length() - K + 1;
    
    // Loop over characters in sequence
    for( int j=0; j < num_substr; j++ ) {
      // Find location in result by looping over K positions
      where = 0;
      for( int k=0; k<K; k++){
        where += seq1[j+k]*Where[k];
      }
      
      // Negative values for alien characters
      if(where >= 0){
        X(where) = true;
      }
    }
    
    // Fill classification vector
    for(int j=0; j < nElem; j++){
      if(X(j)){
        for(int k=0; k < nclass; k++){
          Ci(k) += Q(k,j);          
        }
        X(j) = false;
      }
    }
    
    // Find maximum/maxima
    int ind    = 1;
    double val = Ci(0);
    bool two   = false;
    for(int j=1; j < nclass; j++){
      if( Ci(j) == val ){
        two = true;
      } else {
        if( Ci(j) > val ){
          two = false;
          ind = j+1;
          val = Ci(j);
        }
      }
    }
    
    // Classify
    if(two){
      C(i) = -1;
    } else {
      C(i) = ind;
    }
    
    // Prepare classification vector for next sequence
    for(int j=0; j < nclass; j++){
      Ci(j) = 0;
    }
  }
  
  return C;
}
