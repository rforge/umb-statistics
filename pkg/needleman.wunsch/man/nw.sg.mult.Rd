\name{nw.sg.mult}
\alias{nw.sg.mult}
\alias{nws.sg.mult}
\alias{qnw.sg.mult}
\docType{package}
\title{
Semi-global Needleman-Wunsch for multiple sequenes.
}
\description{
Simple function using Rcpp to compute semi-global Needleman-Wunsch score.
The first set of sequence are short, the second long.
}
\usage{
nw.sg.mult(short, long, S, g)
nws.sg.mult(short, long, S, g)
qnw.sg.mult(short, long)
}
\arguments{
  \item{short}{short sequences (list of integer vectors).}
  \item{long}{long sequences (list of integer vectors).}
  \item{S}{a scoring matrix (\code{nw.sg.mult}) reshaped as vector or a short vector (\code{nws.sg.mult}) based on the scoring matrix (quicker) (see examples).}
  \item{g}{an integer vector holding the gap penalty.}
}
\value{
  \item{D }{Matrix of scores from queries against satellites}
  \item{D_scaled }{Matrix of scores scaled by satellite lengths}
}
\examples{
short <- char2int(c("ACGT", "ACGG"))
long  <- char2int(c("AAACCCACGTATA", "ACGGGACGAGC"))
nw.sg.mult(short, long, matrix(c(1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,1),4,4), -1)
nws.sg.mult(short, long, c(1,-1), -1)
qnw.sg.mult(short, long)
}
