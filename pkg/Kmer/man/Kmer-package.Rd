\name{Kmer-package}
\alias{Kmer-package}
\alias{Kmer}
\alias{Kmer_matrix}
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
system.time(a <- Kmer(c("AAACGTTTT","GTTT"),3,TRUE))
dim(a)
sum(a)
print(a)
system.time(b <- Kmer_classes(c("AAACGTTTT","GTTT"),3,TRUE,c("Fish","Coconut")))
print(b)
}
