# Main K-mer matrix function
Kmer <- function(Seq,K){
  Kmer_matrix(char2int(Seq),K)
}

# Recode character to integer vectors
char2int <- function(str){
  lapply(1:length(str), function(s)as.integer(factor(strsplit(str[s],"")[[1]], levels=c("A","C","G","T")))-1L)
}
