\name{singletons}
\alias{singletons}
\title{
Singletons
}
\description{
Estimates the number of singletons in a rarefied subsample.
}
\usage{
singletons(n,quota)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
  \item{quota}{
An integer specifying the sample size for which the singleton count will be estimated.
}
}
\details{
This function implements a self-evident combinatorial equation that computes the expected number of species found exactly once in a subsample (= singletons), only returning that estimate.
The equation was used by Chao and Jost (2012) in their reformulation of SQS (see \link{sqs}).
}
\references{
Chao, A. and L. Jost. 2012. Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. \emph{Ecology} 93:2533-2547.
}
\author{
John Alroy
}
\examples{
# CEGS distribution
n <- rcegs(100)
singletons(n[n > 0],50)

# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
singletons(n[n > 0],50)

# scaled odds distribution
n <- rodds(100,1)
singletons(n[n > 0],50)
}
