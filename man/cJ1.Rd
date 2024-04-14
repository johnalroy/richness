\name{cJ1}
\alias{cJ1}
\title{
The corrected first-order jackknife richness estimator
}
\description{
Estimates richness using the corrected first-order jackknife equation.
}
\usage{
cJ1(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function implements the corrected first-order jackknife equation of Alroy (2020), only returning the richness estimate.
}
\references{
Alroy, J. 2020. On four measures of taxonomic richness. \emph{Paleobiology} 46:158-175.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
cJ1(n[n > 0])

# scaled odds distribution
n <- rodds(100,1)
cJ1(n[n > 0])

# half-power exponential distribution
n <- rhp(100,1)
cJ1(n[n > 0])
}
