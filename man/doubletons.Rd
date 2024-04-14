\name{doubletons}
\alias{doubletons}
\title{
Doubletons
}
\description{
Estimates the number of doubletons in a rarefied subsample.
}
\usage{
doubletons(n,quota)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
  \item{quota}{
An integer specifying the sample size for which the doubleton count will be estimated.
}
}
\details{
This function implements a self-evident combinatorial equation that computes the expected number of species found exactly twice in a subsample (= doubletons), only returning that estimate.
A similar equation for singletons was used by Chao and Jost (2012) in their reformulation of SQS (see \link{sqs}).
}
\references{
Chao, A. and L. Jost. 2012. Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. \emph{Ecology} 93:2533-2547.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
doubletons(n[n > 0],50)

# scaled odds distribution
n <- rodds(100,1)
doubletons(n[n > 0],50)

# half-power exponential distribution
n <- rhp(100,1)
doubletons(n[n > 0],50)
}
