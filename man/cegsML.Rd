\name{cegsML}
\alias{cegsML}
\title{
The compound exponential-geometric series fit by maximum likelihood
}
\description{
Uses maximum likelihood to fit a species abundance distribution based on a count vector to the compound exponential-geometric series.
}
\usage{
cegsML(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function fits the compound exponential-geometric series (CEGS) distribution, returning the corrected Akaike information criterion (AICc) fit statistic and the rank and species abundance distributions in addition to estimated richness and the distribution's two governing parameters.

CEGS assumes each count is generated by a random draw from the geometric series (GS), but with the governing parameter \emph{p} of the GS being itself a function of a random exponential variate \emph{E} = -log(\emph{U}), where \emph{U} is uniformly distributed between 0 and 1. In other words, each count is doubly random. This is a compound distribution in the standard statistical sense.

The function connects \emph{p} and \emph{E} via two steps. First, because the mean of a GS draw is \emph{p}/(1 - \emph{p}) and the mean of an exponential draw is 1/\eqn{\lambda} where \eqn{\lambda} is a scaling parameter, it follows that \emph{p} = the integral of 1/(\emph{E}/\eqn{\lambda} + 1). Second, we assume that the exponential draw varies according to a power law with a parameter \eqn{\gamma}:

\emph{p} = integral 1/(\emph{E}^\eqn{\gamma}/\eqn{\lambda} + 1)

The computational expression for the probability mass function (PMF) is:

\code{
p <- 1 / ((-log(U))^g / l + 1)
}

\code{
dgeom(x - 1,p)
}

where e = \emph{E}, g = \eqn{\gamma}, l = \eqn{\lambda}, and x = the integer series starting with 1, which corresponds with the observed count classes (e.g., classes 1 and 2 are the respective counts of singletons and doubletons). To obtain a particular probability, the integral is computed over \emph{U}'s range of 0 to 1 by means of the \emph{integrate} function. The equation assumes that all counts are greater than zero (hence x - 1 instead of x).

The model is fit by optimising the AICc using the \emph{stats4} function \emph{mle}, as with some other two-parameter model fitting functions in this package (e.g., \link{weibull}).
The likelihood calculation that yields the AICc is explained in the package details.

\link{cegsB} and \link{cegsLD} produce species richness estimates with lower variance. \link{cegsB} assumes a prior distribution and \link{cegsLD} does not. \emph{cegsML} is faster than either one.

The function \link{rcegs} computes random integers based on the distribution.
}
\value{
  \item{richness}{The sample's estimated species richness.}
  \item{scale}{The estimated \eqn{\lambda} (scale) parameter of the CEGS distribution.}
  \item{shape}{The estimated \eqn{\gamma} (shape) parameter of the CEGS distribution.}
  \item{AICc}{The corrected Akaike information criterion value computed from the log likelihood.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\author{
John Alroy
}
\examples{
# CEGS distribution
n <- rcegs(100)
cegsML(n[n > 0])$richness

# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
cegsML(n[n > 0])$richness

# scaled odds distribution
n <- rodds(100,1)
cegsML(n[n > 0])$richness
}
