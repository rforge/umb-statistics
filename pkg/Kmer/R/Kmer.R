# Main K-mer matrix function
Kmer <- function(Seq, K, names=FALSE){
  X <- Kmer_matrix(char2int(Seq),K,names)
  rownames(X) <- names(Seq)
  X
}

# Main K-mer matrix function with alien characters
KmerAlien <- function(Seq, K, names=FALSE){
  X <- Kmer_matrix_alien(char2intAlien(Seq,K),K,names)
  rownames(X) <- names(Seq)
  X
}

# Main K-mer matrix function (RDP)
KmerRDP <- function(Seq, K, names=FALSE){
  X <- Kmer_matrix_RDP(char2int(Seq),K,names)
  rownames(X) <- names(Seq)
  X
}

# Main K-mer matrix function with alien characters (RDP)
KmerAlienRDP <- function(Seq, K, names=FALSE){
  X <- Kmer_matrix_alien_RDP(char2intAlien(Seq,K),K,names)
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
  }
  X <- Kmer_matrix_classes(char2int(Seq), K, names, classesInt, max(classesInt))
  rownames(X$X) <- names(Seq)
  if(is.character(classes))
    rownames(X$C) <- classLevels
  X
}

# RDP training function
rdpTrain <- function(Seq, classes, K=8, names=FALSE){
  classesInt <- classes
  if(is.character(classes)){
    classesInt  <- factor(classes)
    classLevels <- levels(classesInt)
    classesInt  <- as.integer(classesInt)
  }
  C <- Kmer_matrix_class_alien_RDP(char2intAlien(Seq,K), K, names, classesInt-1, max(classesInt), as.numeric(table(classes)))
  if(is.character(classes))
    rownames(C) <- classLevels
  C
}

# RDP classification
rdpClassify <- function(Seq, Q){
  K <- as.integer( log2( dim( Q )[2] )/2 )
  C <- Kmer_classify_alien_RDP(Q, char2intAlien(Seq,K), K)
  C[C<0] <- NA
  rownames(Q)[C]
}

# Recode character to integer vectors
char2int <- function(str){
  lapply(1:length(str), function(s)as.integer(factor(strsplit(str[s],"")[[1]], levels=c("A","C","G","T")))-1L)
}
char2intAlien <- function(str, K){
  lapply(1:length(str), function(s){
    x <- as.integer(factor(strsplit(str[s],"")[[1]], levels=c("A","C","G","T","U","N","X","W","S","M","K","R","Y","B","D","H","V","-")))-1L
    x[x==4] <- 3
    x[x>4]  <- -4^(K+2)
    x
  })
}
