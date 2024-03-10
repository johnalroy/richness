\name{jackknife2}
\alias{jackknife2}
\title{
The second-order jackknife richness estimator
}
\description{
Estimates richness using the second-order jackknife equation.
}
\usage{
jackknife2(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function implements the second-order jackknife equation of Burnham and Overton (1978), only returning the richness estimate.
}
\references{
Burnham, K. P. and W. S. Overton. 1978. Estimation of the size of a closed population when capture probabilities vary among animals. \emph{Biometrika} 65:625-633.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
jackknife2(n[n > 0])

# scaled odds distribution
n <- rodds(100,1)
jackknife2(n[n > 0])

# exponential-to-e distribution
n <- rexp2e(100,1)
jackknife2(n[n > 0])
}
