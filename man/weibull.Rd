\name{weibull}
\alias{weibull}
\title{
The discrete Weibull distribution
}
\description{
Fits a species abundance distribution based on a count vector to the discrete Weibull distribution.
}
\usage{
weibull(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function fits the discrete Weibull distribution of Nakagawa and Osaki (1975), returning the AICc fit statistic and the rank and species abundance distributions in addition to estimated richness and the distribution's governing parameters (scale and shape). The model is fit by optimising the AICc using the \emph{stats4} function \emph{lme}, as with some of this package's other functions.
The likelihood calculation that yields the AICc is explained in the package details.

This function occasionally crashes because the two-parameter likelihood surface is flat, leading \emph{lme} to consider extreme combinations of values that cause errors.
}
\value{
  \item{richness}{The sample's estimated species richness.}
  \item{scale}{The estimated scale parameter of the Weibull distribution's SAD equation.}
  \item{shape}{The estimated shape parameter of the Weibull distribution's SAD equation.}
  \item{AICc}{The corrected Akaike information criterion value computed from the log likelihood.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\references{
Nakagawa, R., and S. Osaki. 1975. The discrete Weibull distribution. \emph{IEEE Transactions on Reliability} 24: 300–301.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
weibull(n[n > 0])$richness

# scaled odds distribution
n <- rodds(100,1)
weibull(n[n > 0])$richness

# half-power exponential distribution
n <- rhp(100,1)
weibull(n[n > 0])$richness
}
