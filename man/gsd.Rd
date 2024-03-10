\name{gsd}
\alias{gsd}
\title{
The geometric series abundance distribution
}
\description{
Fits a species abundance distribution based on a count vector to the geometric series distribution.
}
\usage{
gsd(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function fits the geometric series distribution to a species abundance distribution (SAD), returning the corrected Akaike information criterion (AICc) fit statistic and the rank and species abundance distributions in addition to the distribution's governing parameter (\emph{k}).
The richness estimate is just the observed number of species divided by 1 minus \emph{k}.
The likelihood calculation that yields the AICc is explained in the package details.

The geometric series is a special one-parameter case of the negative binomial that only allows the probability term to vary.
\link{nbin} allows users to vary the negative binomial's size parameter.

Ecologists are most often familiar with the geometric series based on the work of Motomura (1932), who applied it in a fundamentally different context.
Motomura suggested fitting the series to ranked abundance distributions (RADs).
So if species 1 is the most common and so on, then the predicted value of count \emph{i} is the total number of individuals in the community times (1 - \emph{k})^(\emph{i} - 1) divided by a scaling constant (see He and Tang 2008).

The \emph{gsd} function instead computes the proportion of all species expected to have an exact integer count of 1, 2, 3, etc., i.e., an SAD sensu stricto. All functions in this package that fit distributions also work with SADs. It is very important not to assume that \emph{gsd} directly fits the geometric series to RADs following Motomura's approach.
}
\value{
  \item{richness}{The sample's estimated species richness.}
  \item{k}{The governing parameter \emph{k} of the geometric series distribution.}
  \item{lambda}{An alternate form of the governing parameter, namely \emph{k}/(1 - \emph{k}), that makes sense when the distribution is interpreted as a mixed Poisson exponential.}
  \item{AICc}{The corrected Akaike information criterion value computed from the log likelihood.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\references{
He, F., and D. Tang. 2008. Estimating the niche preemption parameter of the geometric series. \emph{Acta Oecologica} 33:105-107.

Motomura, I. 1932. A statistical treatment of associations. \emph{Japanese Journal of Zoology} 44:379-383.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
gsd(n[n > 0])$richness

# scaled odds distribution
n <- rodds(100,1)
gsd(n[n > 0])$richness

# exponential-to-e distribution
n <- rexp2e(100,1)
gsd(n[n > 0])$richness
}
