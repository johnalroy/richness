\name{pln}
\alias{pln}
\title{
The Poisson log normal abundance distribution
}
\description{
Fits a species abundance distribution based on a count vector to the Poisson log normal distribution.
}
\usage{
pln(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function returns statistics based on the Poisson log normal distribution described by Bulmer (1974), returning the corrected Akaike information criterion (AICc) fit statistic and the rank and species abundance distributions in addition to estimated richness and the distribution's governing parameters \emph{mu} and \emph{sigma}.
It depends on the function \emph{poilogMLE} in the \emph{sads} library to compute the parameters.
}
\value{
  \item{richness}{The sample's estimated species richness.}
  \item{mu}{The underlying mean of the fitted log normal distribution.}
  \item{sigma}{The underlying standard deviation of the fitted log normal distribution.}
  \item{AICc}{The corrected Akaike information criterion value computed from the log likelihood.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\references{
Bulmer, M. G. 1974. On fitting the Poisson lognormal distribution to species abundance data. \emph{Biometrics} 30:651-660.
}
\author{
John Alroy
}
\examples{
# CEGS distribution
n <- rcegs(100)
pln(n[n > 0])$richness

# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
pln(n[n > 0])$richness

# scaled odds distribution
n <- rodds(100,1)
pln(n[n > 0])$richness
}
