\name{Anova.lm}
\alias{Anova.lm}
\alias{Anova}
\encoding{latin1}
\title{
Analysis of variance with SS type II or III (including mixed models).
}
\description{
Replacement for \code{Anova.lm} in package \code{car}.
}
\usage{
\method{Anova}{lm}(mod, type="III", ...)
}
\arguments{
  \item{mod}{fitted model, e.g. from \code{lm} or \code{lmer}.}
  \item{type}{SS type indicated by 'II' (default), 'III', 2 or 3.}
  \item{\dots}{possible additional arguments passed to the underlying functions.}
}
\value{
	Returns appropriate analysis of variance or halts if unsupported input is detected.
}
\author{
Kristian Hovde Liland
}

\seealso{
\code{\link{print.AnovaMix}}, \code{\link{AnovaMix}}, \code{\link{lm}}
}
\examples{
mixlm <- lm(y~x*r(z),
		    data = data.frame(y = rnorm(8),
							  x = factor(c(rep(1,4),rep(0,4))),
							  z = factor(rep(c(1,0),4))))
Anova(mixlm, type="III")
}
\keyword{ ANOVA }
