\name{shannon}
\alias{shannon}
\title{
Shannon's H
}
\description{
Computes Shannon's \emph{H}.
}
\usage{
shannon(n,exponentiate=TRUE)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
  \item{exponentiate}{
Whether to exponentiate the returned value.
}
}
\details{
This function implements the equation of Shannon (1948), only returning the \emph{H} estimate.
The exponential of \emph{H} is a Hill number (Hill 1973) that converges on true species richness as an abundance distribution approaches uniformity.
}
\references{
Hill, M. O. 1973. Diversity and evenness: a unifying notation and its consequences. \emph{Ecology} 54:427-432.

Shannon, C. E. 1948. A mathematical theory of communication. \emph{Bell System Technical Journal} 27:379-423.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
shannon(n[n > 0])

# scaled odds distribution
n <- rodds(100,1)
shannon(n[n > 0])

# exponential-to-e distribution
n <- rexp2e(100,1)
shannon(n[n > 0])
}
