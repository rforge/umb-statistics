\name{char2int}
\alias{char2int}
\docType{package}
\title{
Character to integer.
}
\description{
Convert character coded nucleotide sequences to a list of integer vectors.
A, C, G, and T are coded as 0, 1, 2, and 3, respectively.
}
\usage{
char2int(str)
}
\arguments{
  \item{str}{a character vector of sequences.}
}
\examples{
char2int("ACTGCTAG")
}
