\name{gsindex}
\alias{gsindex}
\title{
The geometric series index richness estimator
}
\description{
Estimates richness using the geometric series index.
}
\usage{
gsindex(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function implements the geometric series index of Kerr and Alroy (2023), only returning the richness estimate.
}
\references{
Kerr, M. R. and J. Alroy. 2023. Body size and abundance are decoupled from species richness in Australian marine bivalves. \emph{Frontiers of Biogeography} 15(4). http://dx.doi.org/10.21425/F5FBG58651
}
\author{
John Alroy
}
\examples{
# CEGS distribution
n <- rcegs(100)
gsindex(n[n > 0])

# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
gsindex(n[n > 0])

# scaled odds distribution
n <- rodds(100,1)
gsindex(n[n > 0])
}
