\name{bsd}
\alias{bsd}
\title{
The broken stick abundance distribution
}
\description{
Fits a species abundance distribution based on a count vector to the broken stick distribution.
}
\usage{
bsd(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function fits the distribution described by MacArthur (1957), returning the corrected Akaike information criterion (AICc) fit statistic and the rank and species abundance distributions in addition to a richness estimate.
The estimate determines the other returned values. 
As originally described, observed richness was assumed to be a fixed parameter and not extrapolated.
Here, richness is estimated by optimising the AICc using the \emph{stats4} function \emph{lme}, as with this package's other functions.
The likelihood calculation that yields the AICc is explained in the package details.

The broken stick model is so poor that it may not be possible to fit it sensiblyby likelihood. Specifically, a returned richness value of two times raw richness most often suggests that no good solution was found.
}
\value{
  \item{richness}{The sample's estimated species richness.}
  \item{AICc}{The corrected Akaike information criterion value computed from the log likelihood.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\references{
MacArthur, R. H. 1957. On the relative abundance of bird species. \emph{Proceedings of the National Academy of Sciences USA} 43:293-295.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
bsd(n[n > 0])$richness

# scaled odds distribution
n <- rodds(100,1)
bsd(n[n > 0])$richness

# half-power exponential distribution
n <- rhp(100,1)
bsd(n[n > 0])$richness
}
