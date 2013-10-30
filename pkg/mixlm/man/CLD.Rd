\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{cld}
\alias{cld}
\alias{print.cld}
\title{
	Compact letter display for pairwise multiple comparisons.
}
\description{
cld takes objects of classes TukeyHSD, Fisher and Bonferroni.
}
\usage{
cld(object, level = 1 - attr(object, "conf.level"))
}
\arguments{
  \item{object}{
	\code{object} of class TukeyHSD, Fisher or Bonferroni.
}
  \item{level}{
	\code{level} of pairwise tests.
}
}
\details{
	Basic compact letter display. \code{print} method
	\code{print.cld} also included.
}
\value{
  \item{comp }{logical matrix of group memberships.}
  \item{cld }{character matrix of letters.}
}
\author{
Kristian Hovde Liland
}

\seealso{
\code{\link{help}}, ~~~
  \code{\link[multcomp]{cld}} in package \CRANpkg{multcomp}.
}
\examples{
summary(fm1 <- lm(breaks ~ wool + tension, data = warpbreaks))
TukeyHSD(fm1, "tension", ordered = TRUE)
cld(TukeyHSD(fm1, "tension", ordered = TRUE))
}
