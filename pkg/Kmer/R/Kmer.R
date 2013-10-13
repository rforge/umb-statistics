# Main K-mer matrix function
Kmer <- function(Seq, K, names=FALSE){
  X <- Kmer_matrix(char2int(Seq),K,names)
  rownames(X) <- names(Seq)
  X
}

# Recode character to integer vectors
char2int <- function(str){
  lapply(1:length(str), function(s)as.integer(factor(strsplit(str[s],"")[[1]], levels=c("A","C","G","T")))-1L)
}
