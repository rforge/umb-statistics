\name{Kmer-package}
\alias{Kmer-package}
\alias{Kmer}
\alias{KmerRDP}
\alias{KmerAlien}
\alias{KmerAlienRDP}
\alias{Kmer_classes}
\alias{Kmer_matrix}
\alias{Kmer_matrix_RDP}
\alias{Kmer_matrix_alien}
\alias{Kmer_matrix_alien_RDP}
\alias{Kmer_matrix_classes}
\alias{char2int}
\docType{package}
\title{
  K-mer calculation using C++ with Rcpp.
}
\description{
  Calculates K-mers from string vectors.
}
\details{
  No details so far.
}
\author{
  Kristian Hovde Liland
}
\references{
  TODO:
}
\seealso{
  Please don't.
}
\examples{
# Basic
system.time(a <- Kmer(c("AAACGTTTT","GTTT"),3,TRUE))
dim(a)
sum(a)
print(a)

# 0/1 instead of count
system.time(aR <- KmerRDP(c("AAACGTTTT","GTTT"),3,TRUE))

# Alien characters
system.time(aA <- KmerAlien(c("AAANGTTTT","GTTT"),3,TRUE))
system.time(aAR <- KmerAlienRDP(c("AAANGTTTT","GTTT"),3,TRUE))

# Count for classes
system.time(b <- Kmer_classes(c("AAACGTTTT","GTTT"),3,TRUE,c("Fish","Coconut")))
print(b)
}
