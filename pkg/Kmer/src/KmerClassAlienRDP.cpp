#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Kmer_matrix_class_alien_RDP( SEXP seqs, int K, bool names, SEXP classesIn, int nclass, SEXP Min ) {
  
  Rcpp::List strings(seqs);
  IntegerVector classes(classesIn);
  NumericVector M(Min);
  int where = 0;
  int nElem = pow(4,K);
  int num_strings = strings.length();
  double pj = 0;
  NumericMatrix C(nclass, nElem);
  NumericVector p(nElem);
  LogicalVector X(nElem);
  std::vector<int> Where(K);
  
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
    
    // Summarize for groups, e.g. genera
    for(int j=0; j < nElem; j++){
      if(X(j)){
        ++C(classes[i],j);
        ++p(j);
        X(j) = false;
      }
    }
  }
  
  // Convert from counts to pseudo probabilities
  for(int j=0; j < nElem; j++){
    pj = (p(j)+0.5)/(num_strings+1);
    for(int i=0; i < nclass; i++){
      C(i,j) = (C(i,j)+pj)/(M(i)+1.0);
    }
  }  
  
  // Create dimnames for output matrix
  if(names){
    int N = pow(4,K);
    Rcpp::CharacterVector ACGT = Rcpp::CharacterVector::create("A","C","G","T");
    Rcpp::CharacterVector ACGTs(N);
    std::vector< std::vector< std::string > > matr;
    matr.resize( K , std::vector<std::string>( N ) );
    Rcpp::CharacterVector cnms(N);
    for(int i=0; i<K; i++){
      ACGTs = rep(rep_each(ACGT, pow(4,(K-1-i))), pow(4,i));
      matr[i] = Rcpp::as< std::vector< std::string > >(ACGTs);
    }
    for(int i=0; i<N; i++){
      std::stringstream ss;
      for(int j=0; j<K; j++){
        ss << matr[j][i];
      }
      cnms[i] = ss.str();      
    }
    Rcpp::List dimnms = Rcpp::List::create(R_NilValue, cnms);
    C.attr("dimnames") = dimnms;
  }
  
  return C;
}
