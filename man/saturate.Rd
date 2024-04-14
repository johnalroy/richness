\name{saturate}
\alias{saturate}
\title{
Saturated abundance distributions
}
\description{
Uses bootstrapping to compute a probability mass function based on a supplied list of counts and on the assumption of a geometric sampling process.
}
\usage{
saturate(n,trials=100,dgeoms=NA)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
  \item{quota}{
An integer specifying the number of bootstrapping trials.
}
  \item{dgeoms}{
The matrix of probabilities returned by \link{dgeoms}.
}
}
\details{
This function fits a saturated model to count data. It returns a smoothed list of probabilities that specify a species abundance distribution (SAD).

The goal is to test whether theoretical abundance distributions can be rejected because they do not fit the data as well as a saturated model.
This can be done by comparing the log likelihoods produced by this function with those produced by the fitted distributions.
See the example below.

Smoothing is accomplished by (1) bootstrapping the list of counts, and (2) assuming that each bootstrapped count is generated separately by a geometric sampling process. The geometric series density vectors matched to the counts are averaged to produce the SAD.
Bootstrapping is carried out by sampling the count vector with replacement.
}
\value{
  \item{LL}{The mean log likelihood across all of the bootstrap trials.}
  \item{LLs}{A sorted list of the log likelihoods produced during the bootstrap trials.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\author{
John Alroy
}
\examples{
# prepare the probability matrix
dg <- dgeoms()

# generate a random list of non-zero counts from the scaled odds distribution
n <- rodds(100)
n <- n[n > 0]

# fit the bootstrapped distribution
s <- saturate(n,100,dg)

# plot the data
plot(rev(sort(n)),type='l',log='y',xlab='Rank of count',ylab='Count of individuals')
lines(rev(s$fitted.RAD),col='red')

# compare the log likelihoods to the observed values for fitted model
s$mean.LL
AICc2LL(n,sodds(n)$AICc)
AICc2LL(n,halfpower(n)$AICc)
AICc2LL(n,logd(n)$AICc)
}
