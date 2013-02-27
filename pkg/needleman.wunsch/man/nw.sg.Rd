\name{nw.sg}
\alias{nw.sg}
\alias{nws.sg}
\alias{qnw.sg}
\docType{package}
\title{
Semi-global Needleman-Wunsch.
}
\description{
Simple function using Rcpp to compute semi-global Needleman-Wunsch score.
The first sequence is short, the second long.
}
\usage{
nw.sg(short, long, S, g)
nws.sg(short, long, S, g)
qnw.sg(short, long)
}
\arguments{
  \item{short}{a short integer vector coding "ACGT" as c(0,1,2,3).}
  \item{long}{a long integer vector coding "ACGT" as c(0,1,2,3).}
  \item{S}{a scoring matrix (\code{nw.sq}) reshaped as vector or a short vector (\code{nws.sq}) based on the scoring matrix (quicker) (see examples).}
  \item{g}{an integer vector holding the gap penalty.}
}
\value{
  \item{D}{Score value of pairwise comparison.}
  \item{L}{Length of pairwise alignment.}
  \item{j}{Endpoint of short sequence in long sequence.}
}
\examples{
short <- char2int("ACTCCTCTAAGAG")
long <- char2int("ACTATATCGGCGCACCCTATAAGAGCACGCGAGA")
nw.sg(short[[1]], long[[1]], c(1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,1), -1)
nws.sg(short[[1]], long[[1]], c(1,-1), -1)
qnw.sg(short[[1]], long[[1]])
}
