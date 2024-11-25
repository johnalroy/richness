\name{chao1}
\alias{chao1}
\title{
The Chao 1 richness estimator
}
\description{
Estimates richness using the Chao 1 equation.
}
\usage{
chao1(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function implements the Chao 1 equation of Chao (1984, 1987), only returning the richness estimate.
}
\references{
Chao, A. 1984. Nonparametric estimation of the number of classes in a population. \emph{Scandinavian Journal of Statistics} 11, 265-270.

Chao, A. 1987. Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics} 43:783-791.
}
\author{
John Alroy
}
\examples{
# CEGS distribution
n <- rcegs(100)
chao1(n[n > 0])

# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
chao1(n[n > 0])

# scaled odds distribution
n <- rodds(100,1)
chao1(n[n > 0])
}
