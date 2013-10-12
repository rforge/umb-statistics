#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix Kmer_matrix( SEXP seqs, int K ) {
  
  Rcpp::List strings(seqs);
  int num_strings = strings.length();
  IntegerMatrix X(num_strings, pow(4,K));
  int where = 0;
  std::vector<int> Where(K);
  
  for(int i=0; i < K; i++){
    Where[i] = pow(4,K-i-1);
  }
  
  for( int i=0; i < num_strings; i++ ) {
    SEXP s1 = strings[i];
    Rcpp::IntegerVector seq1(s1);
    int num_substr = seq1.length() - K + 1;
    
    for( int j=0; j < num_substr; j++ ) {
      where = 0;
      for( int k=0; k<K; k++){
        where += seq1[j+k]*Where[k];
      }
      
      ++X(i, where);
    }
  }
  
  return X;
}
