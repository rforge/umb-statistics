\name{nw}
\alias{nw}
\alias{nws}
\alias{qnw}
\docType{package}
\title{
Needleman-Wunsch.
}
\description{
Simple function using Rcpp to compute Needleman-Wunsch score.
}
\usage{
nw(seq1, seq2, S, g)
nws(seq1, seq2, S, g)
qnw(seq1, seq2)
}
\arguments{
  \item{seq1}{an integer vector coding "ACGT" as c(0,1,2,3).}
  \item{seq2}{an integer vector coding "ACGT" as c(0,1,2,3).}
  \item{S}{a scoring matrix (\code{nw}) reshaped as vector or a short vector (\code{nws}) based on the scoring matrix (quicker) (see examples).}
  \item{g}{an integer vector holding the gap penalty.}
}
\value{
  \item{D}{Score value of pairwise comparison.}
  \item{L}{Length of pairwise alignment.}
}
\examples{
seq1 <- char2int("ACTCCTCTAAGAG")
seq2 <- char2int("ACCCTATAAGAGC")
nw(seq1[[1]], seq2[[1]], c(1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,1), -1)
nws(seq1[[1]], seq2[[1]], c(1,-1), -1)
qnw(seq1[[1]], seq2[[1]])
}
