\name{zsm}
\alias{zsm}
\title{
The zero-sum multinomial abundance distribution
}
\description{
Fits a species abundance distribution based on a count vector to the zero-sum multinomial distribution.
}
\usage{
zsm(n)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
}
\details{
This function returns statistics based on the zero-sum multinomial distribution described by Hubbell (1997, 2001), including the corrected Akaike information criterion (AICc) fit statistic and the rank and species abundance distributions plus the distribution's two governing parameters \emph{J} and \emph{theta}.
It uses the function \emph{fitmzsm} in the library \emph{sads} to obtain the two parameters.
}
\value{
  \item{J}{The fitted \emph{J} governing parameter of the distribution.}
  \item{theta}{The fitted \emph{theta} governing parameter of the distribution, which is often close to the log series parameter alpha (see \link{fisher} and \link{logd}.)}
  \item{AICc}{The corrected Akaike information criterion value computed from the log likelihood.}
  \item{fitted.RAD}{The predicted rank abundance distribution.}
  \item{fitted.SAD}{The predicted species abundance distribution.}
}
\references{
Hubbell, S. P. 1997. A unified theory of biogeography and relative species abundance and its application to tropical rain forests and coral reefs. \emph{Coral Reefs} 16:S9-S21.

Hubbell, S. P. 2001. The unified neutral theory of biodiversity and biogeography. Princeton Univ. Press, Princeton, NJ.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
zsm(n[n > 0])$fitted.RAD

# scaled odds distribution
n <- rodds(100,1)
zsm(n[n > 0])$fitted.RAD

# exponential-to-e distribution
n <- rexp2e(100,1)
zsm(n[n > 0])$fitted.RAD
}
