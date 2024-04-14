\name{sqs}
\alias{sqs}
\title{
Shareholder quorum subsampling
}
\description{
Estimates subsampled richness using shareholder quorum subsampling (coverage-based rarefaction).
}
\usage{
sqs(n,quorum)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
  \item{quorum}{
A desired (target) level of frequency distribution coverage sensu Good (1953).
}
}
\details{
This function implements the shareholder quorum subsampling method of Alroy (2010).
It uses the analytical solution described by Chao and Jost (2012), who called the method coverage-based rarefaction.
Output statistics are computed by finding the number of individuals that must be drawn from the sample to reach the target coverage level (= quorum).
}
\value{
  \item{raw.richness}{The number of species in the full sample.}
  \item{individuals}{The number of individuals in the full sample.}
  \item{subsampled.richness}{The number of species expected to be found at the quorum.}
  \item{subsampled.singletons}{The number of singletons expected to be found at the quorum.}
  \item{individuals.drawn}{The number of individuals that need to be drawn on average to reach the quorum.}
}
\references{
Alroy, J. 2010. The shifting balance of diversity among major animal groups. \emph{Science} 329:1191-1194.

Chao, A. and L. Jost. 2012. Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. \emph{Ecology} 93:2533-2547.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
sqs(n[n > 0],0.5)

# scaled odds distribution
n <- rodds(100,1)
sqs(n[n > 0],0.5)

# half-power exponential distribution
n <- rhp(100,1)
sqs(n[n > 0],0.5)
}
