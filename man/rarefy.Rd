\name{rarefy}
\alias{rarefy}
\title{
Rarefaction
}
\description{
Estimates subsampled richness using rarefaction.
}
\usage{
rarefy(n,quota)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
  \item{quota}{
An integer specifying the sample size for which the subsampled richness will be estimated.
}
}
\details{
This function implements the analytical rarefaction equation of Hurlbert (1971), which was inspired by the rarefaction algorithm of Sanders (1968). It only returns the richness estimate.
}
\references{
Sanders, H. L. 1968. Marine benthic diversity: a comparative study. \emph{American Naturalist} 102:243-282.

Hurlbert, S. H. 1971. The nonconcept of species diversity: a critique and alternative parameters. \emph{Ecology} 52:577-586.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
rarefy(n[n > 0],50)

# scaled odds distribution
n <- rodds(100,1)
rarefy(n[n > 0],50)

# exponential-to-e distribution
n <- rexp2e(100,1)
rarefy(n[n > 0],50)
}
