\name{fisher}
\alias{fisher}
\title{
Fisher's alpha
}
\description{
Estimates Fisher's alpha.
}
\usage{
fisher(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function computes the log series parameter alpha described by Fisher et al. (1943) using the recursive equation implied by Chatfield (1969) and detailed by May (1975).
The \emph{sads} library uses the same equation.
The function only returns the estimate.
}
\references{
Chatfield, C. 1969. On estimating the parameters of the logarithmic series and negative binomial distributions. \emph{Biometrika} 56:411-414.

Fisher, R. A., A. S. Corbet, and C. B. Williams. 1943. The relation between the number of species and the number of individuals in a random sample of an animal population. \emph{Journal of Animal Ecology} 12:42-58.

May, R. M. 1975. Patterns of species abundance and diversity. Pp. 81-120 in M. L. Cody and J. M. Diamond (rds.), \emph{Ecology and Evolution of Communities}. Belknap Press.
}
\author{
John Alroy
}
\examples{
# CEGS distribution
n <- rcegs(100)
fisher(n[n > 0])

# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
fisher(n[n > 0])

# scaled odds distribution
n <- rodds(100,1)
fisher(n[n > 0])
}
