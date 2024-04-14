\name{nbin}
\alias{nbin}
\title{
The negative binomial distribution
}
\description{
Fits a species abundance distribution based on a count vector to the negative binomial distribution.
}
\usage{
nbin(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function fits the negative binomial distribution, returning the AICc fit statistic and the rank and species abundance distributions in addition to estimated richness and the distribution's governing parameters (size and probability). The model is fit by optimising the AICc using the \emph{stats4} function \emph{lme}, as with this package's other functions.
The likelihood calculation that yields the AICc is explained in the package details.

The geometric series is a special case of the negative binomial in which the size parameter is fixed at 1. The \link{gsd} function allows fitting this version of the model.
}
\value{
  \item{richness}{The sample's estimated species richness.}
  \item{size}{The estimated size parameter of the negative binomial distribution's SAD equation.}
  \item{probability}{The estimated probability parameter of the negative binomial distribution's SAD equation.}
  \item{AICc}{The corrected Akaike information criterion value computed from the log likelihood.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
nbin(n[n > 0])$fitted.RAD

# scaled odds distribution
n <- rodds(100,1)
nbin(n[n > 0])$fitted.RAD

# half-power exponential distribution
n <- rhp(100,1)
nbin(n[n > 0])$fitted.RAD
}
