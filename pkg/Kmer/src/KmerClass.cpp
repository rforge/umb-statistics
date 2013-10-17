#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List Kmer_matrix_classes( SEXP seqs, int K, bool names, SEXP classesIn, int nclass ) {
  
  Rcpp::List strings(seqs);
  Rcpp::IntegerVector classes(classesIn);
  int num_strings = strings.length();
  IntegerMatrix X(num_strings, pow(4,K));
  IntegerMatrix C(nclass, pow(4,K));
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
      ++C(classes[i]-1,where);
    }
  }
  
  
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
    X.attr("dimnames") = dimnms;
    C.attr("dimnames") = dimnms;
  }
  
  return Rcpp::List::create(Rcpp::Named("X") = X, Rcpp::Named("C") = C);
}
