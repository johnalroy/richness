\name{logd}
\alias{logd}
\title{
The log series abundance distribution
}
\description{
Fits a species abundance distribution based on a count vector to the log series distribution.
}
\usage{
logd(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function fits the log series distribution described by Fisher et al. (1943), returning the corrected Akaike information criterion (AICc) fit statistic and the rank and species abundance distributions in addition to the distribution's governing parameter (Fisher's alpha).
A meaningful richness estimate cannot be obtained because the distribution does not indicate what proportion of species is unsampled.
The function relies upon \link{fisher} to compute alpha using the standard recursive algorithm implied by Chatfield (1969) and detailed by May (1975).
}
\value{
  \item{alpha}{The governing parameter of the log series distribution.}
  \item{AICc}{The corrected Akaike information criterion value computed from the log likelihood.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\references{
Chatfield, C. 1969. On estimating the parameters of the logarithmic series and negative binomial distributions. \emph{Biometrika} 56:411-414.

Fisher, R. A., A. S. Corbet, and C. B. Williams. 1943. The relation between the number of species and the number of individuals in a random sample of an animal population. \emph{Journal of Animal Ecology} 12:42-58.

May, R. M. 1975. Patterns of species abundance and diversity. Pp. 81-120 in M. L. Cody and J. M. Diamond (rds.), \emph{Ecology and Evolution of Communities}. Belknap Press.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
logd(n[n > 0])$fitted.RAD

# scaled odds distribution
n <- rodds(100,1)
logd(n[n > 0])$fitted.RAD

# half-power exponential distribution
n <- rhp(100,1)
logd(n[n > 0])$fitted.RAD
}
