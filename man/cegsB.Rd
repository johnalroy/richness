\name{cegsB}
\alias{cegsB}
\title{
The compound exponential-geometric series fit by a Bayesian calculation
}
\description{
Uses a Bayesian calculation to fit a species abundance distribution to the compound exponential-geometric series.
}
\usage{
cegsB(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function fits the compound exponential-geometric series (CEGS) distribution, returning the corrected Akaike information criterion (AICc) fit statistic and the rank and species abundance distributions in addition to estimated richness and the distribution's two governing parameters.

CEGS assumes each count is generated by a random draw from the geometric series (GS), but with the governing parameter \emph{p} of the GS being itself a function of a random exponential variate \emph{E} = -log(\emph{U}), where \emph{U} is uniformly distributed between 0 and 1. In other words, each count is doubly random. The expression is fitted to data by varying a scale term \eqn{\lambda} and a shape term \eqn{\gamma}, which feature in the probability mass function at different points.

Further details that include the derivation of the compound distribution are given on the \link{cegsML} manual page.

The Bayesian calculation of \eqn{\lambda} assumes that the priori distribution is exponential. In other words, it is produced by taking negative logs of fractions that are evenly spaced between 0 and 1.
Fifty fractions are used, so there are limits of -log(50/51) = 0.020 to -log(1/51) = 3.932.

For \eqn{\gamma}, the prior distribution is the log ratio of an odds distribution: -log(\emph{U} / (1 - \emph{U})). The result is a symmetrical distribution with moderate variance that is centred on zero and has analogous extreme limits: -3.932 and 3.932.

A grid computation is used to obtain the likelihoods of the count distributions given possible combinations of the parameters (i.e., the conditional probabilities). To do this, a matrix of \eqn{\lambda} quantiles spaced equally from 0 to 1 along the columns and fixed across the rows is created. An orthogonal matrix for \eqn{\gamma} is also created. Likelihoods are found for each combination (= grid point). The matrix of likelihoods is scaled so it sums to 1; it is multiplied across the \eqn{\lambda} quantile matrix; the sum is taken; and finally, the negative log of the sum is taken to obtain the posterior estimate of \eqn{\lambda}. \eqn{\gamma} is calculated in the same way, except that the log of an odds distribution is taken.

An AICc is also reported for comparability with other functions. \link{cegsML} will produce lower AICcs if a likelihood approach is preferred. The likelihood calculation that yields the AICc is explained in the package details.

Like \link{cegsLD}, \emph{cegsB} produces species richness estimates with low variance because it takes all likelihoods into account. However, \link{cegsLD} does not assume a prior distribution. \link{cegsML} is faster than either one.

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
cegsB(n[n > 0])$richness

# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
cegsB(n[n > 0])$richness

# scaled odds distribution
n <- rodds(100,1)
cegsB(n[n > 0])$richness
}
