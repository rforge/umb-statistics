\name{nw.sat.quer}
\alias{nw.sat.quer}
\alias{nws.sat.quer}
\alias{qnw.sat.quer}
\docType{package}
\title{
Needleman-Wunsch for many sequenes.
}
\description{
Few/many to many pairwise comparisons with two sets of sequences
(satellites and queries)
}
\usage{
nw.sat.quer(satellites, queries, S, g)
nws.sat.quer(satellites, queries, S, g)
qnw.sat.quer(satellites, queries)
}
\arguments{
  \item{satellites}{satellite sequences (list of integer vectors).}
  \item{queries}{other sequences (list of integer vectors).}
  \item{S}{a scoring matrix (\code{nw.sat.quer}) reshaped as vector or a short vector (\code{nws.sat.quer}) based on the scoring matrix (quicker) (see examples).}
  \item{g}{an integer vector holding the gap penalty.}
}
\value{
  \item{D }{Matrix of scores from queries against satellites}
  \item{D_scaled }{Matrix of scores scaled by satellite lengths}
}
\examples{
satellites <- char2int(c("ACGT", "ACGG"))
queries    <- char2int(c("ACGT", "ACGGG"))
nw.sat.quer(satellites, queries, c(1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,1), -1)
nws.sat.quer(satellites, queries, c(1,-1), -1)
qnw.sat.quer(satellites, queries)
}
