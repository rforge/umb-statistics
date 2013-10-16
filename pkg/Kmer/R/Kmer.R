# Main K-mer matrix function
Kmer <- function(Seq, K, names=FALSE){
  X <- Kmer_matrix(char2int(Seq),K,names)
  rownames(X) <- names(Seq)
  X
}

# K-mer matrix and sum per class, e.g. per genus
Kmer_classes <- function(Seq, K, names=FALSE, classes){
  classesInt <- classes
  if(is.character(classes)){
    classesInt  <- factor(classes)
    classLevels <- levels(classesInt)
    classesInt  <- as.integer(classesInt)
    classUnique <- unique(classesInt)
  }
  X <- Kmer_matrix_classes(char2int(Seq), K, names, classesInt, max(classesInt))
  rownames(X$X) <- names(Seq)
  if(is.character(classes))
    rownames(X$C) <- classLevels
  X
}

# Recode character to integer vectors
char2int <- function(str){
  lapply(1:length(str), function(s)as.integer(factor(strsplit(str[s],"")[[1]], levels=c("A","C","G","T")))-1L)
}
