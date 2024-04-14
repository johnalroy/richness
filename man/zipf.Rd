\name{zipf}
\alias{zipf}
\title{
The Zipf abundance distribution
}
\description{
Fits a species abundance distribution based on a count vector to the Zipf distribution.
}
\usage{
zipf(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function returns statistics based on the Zipf distribution (see Newman 2005), returning the AICc fit statistic and the rank and species abundance distributions in addition to the distribution's governing parameter (the power).
A meaningful richness estimate cannot be obtained because the distribution assumes richness is infinite.
The function uses uses \emph{fitpower} in the \emph{sads} library to obtain the power.
}
\value{
  \item{exponent}{The fitted governing parameter of the distribution.}
  \item{AICc}{The corrected Akaike information criterion value computed from the log likelihood.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\references{
Newman, M. E. J. 2005. Power laws, Pareto distributions and Zipf's law. \emph{Contemporary Physics} 46:323-351.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
zipf(n[n > 0])$fitted.RAD

# scaled odds distribution
n <- rodds(100,1)
zipf(n[n > 0])$fitted.RAD

# half-power exponential distribution
n <- rhp(100,1)
zipf(n[n > 0])$fitted.RAD
}
