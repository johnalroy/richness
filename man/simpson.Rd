\name{simpson}
\alias{simpson}
\title{
Simpson's D
}
\description{
Estimates richness using the equation.
}
\usage{
simpson(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function implements the corrected Simpson's \emph{D} equation of Hurlbert (1971), only returning the inverse of the \emph{D} estimate.
The original equation was published by Simpson (1949).
The inverse of \emph{D} is a Hill number (Hill 1973) that converges on true species richness as an abundance distribution approaches uniformity.
}
\references{
Hill, M. O. 1973. Diversity and evenness: a unifying notation and its consequences. \emph{Ecology} 54:427-432.

Hurlbert, S. H. 1971. The nonconcept of species diversity: a critique and alternative parameters. \emph{Ecology} 52:577-586.

Simpson, E. H. 1949. Measurement of diversity. \emph{Nature} 163:688.
}
\author{
John Alroy
}
\examples{
# CEGS distribution
n <- rcegs(100)
simpson(n[n > 0])

# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
simpson(n[n > 0])

# scaled odds distribution
n <- rodds(100,1)
simpson(n[n > 0])
}
