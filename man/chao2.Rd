\name{chao2}
\alias{chao2}
\title{
The Chao 2 richness estimator
}
\description{
Estimates richness using the Chao 2 equation.
}
\usage{
chao2(m,s)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of presences per species across a set of samples.
}
  \item{s}{
The number of samples that contributed to the count vector.
}
}
\details{
This function implements the Chao 2 equation of Chao (1987), only returning the richness estimate.
}
\references{
Chao, A. 1987. Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics} 43:783-791.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
m <- matrix(rpois(1000,exp(rnorm(100,sd=2))),100,10)
m <- m[rowSums(m) > 0,]
m[m > 1] <- 1
chao2(rowSums(m),ncol(m))

# scaled odds distribution
n <- rodds(100,1)
m <- matrix(NA,100,10)
for (i in 1:100)
	m[i,] <- rmultinom(1:10,n[i],rep(0.1,10))
m <- m[rowSums(m) > 0,]
m[m > 1] <- 1
chao2(rowSums(m),ncol(m))

# exponential-to-e distribution
m <- matrix(rexp2e(1000,1),100,10)
m <- m[rowSums(m) > 0,]
m[m > 1] <- 1
chao2(rowSums(m),ncol(m))
}
