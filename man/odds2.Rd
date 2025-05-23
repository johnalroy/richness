\name{odds2}
\alias{odds2}
\title{
The odds-2 distribution fit by maximum likelihood
}
\description{
Uses maximum likelihood to fit a species abundance distribution based on a count vector to the odds-2 distribution.
}
\usage{
odds2(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function fits the odds-2 distribution, returning the corrected Akaike information criterion (AICc) fit statistic and the rank and species abundance distributions in addition to estimated richness and the distribution's single governing parameter.

odds-2 is a variant of the scaled odds distribution (see \link{sodds}). It assumes each count is generated by rounding down a value \emph{x} drawn from random odds that are centered with a constant \emph{mu} (m in the code). Here, the odds are based on a square root-transformed random uniform variate \emph{U}. The equation underpinning the distribution is \emph{x} = \emph{mu} (1 - \emph{U}^0.5)/\emph{U}^0.5. The probability mass function follows from rearranging this expression algebraically to solve for \emph{U}. Unlike the scaled odds, the distribution often generates medium-scale variation in counts.

As with the scaled odds (\link{sodds}) and half-power exponential (\link{halfpower}) distributions, the odds-2 distribution can be matched precisely with data derived from a simple birth-death simulation. These models are fully neutral except that they assume individual species have birth and death rates that are fixed and correlated. They are expected to be realistic whenever traits determine demographic properties. 

The analytical expression for the probability mass function (PMF) is lengthy. The algorithmic expression is:

\code{
-diff(1/(m + x)^2) * (m + 1)^2
}

where m = \emph{mu} and x = the positive integer series starting with 1, which corresponds with the count classes (e.g., classes 1 and 2 are the respective counts of singletons and doubletons). The equation assumes that all counts are greater than zero.

The model is fit by optimising the AICc using \emph{optimise}, as with most of this package's other functions.
The likelihood calculation that yields the AICc is explained in the package details.

The function \link{rodds2} computes random integers based on the distribution.
}
\value{
  \item{richness}{The sample's estimated species richness.}
  \item{mu}{The estimated governing parameter of the odds-2 distribution.}
  \item{AICc}{The corrected Akaike information criterion value computed from the log likelihood.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\author{
John Alroy
}
\examples{
# CEGS distribution
n <- rcegs(100)
odds2(n[n > 0])$richness

# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
odds2(n[n > 0])$richness

# scaled odds distribution
n <- rodds(100,1)
odds2(n[n > 0])$richness

# odds-2 distribution
n <- rodds2(100,3)
odds2(n[n > 0])$richness
}
