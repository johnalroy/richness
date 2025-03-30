\name{rcegs}
\alias{rcegs}
\title{
Random generation from the CEGS distribution
}
\description{
Draws counts at random from the CEGS distribution.
}
\usage{
rcegs(n=100,l=1,g=2)
}
\arguments{
  \item{n}{
The number of counts to be drawn.
}
  \item{l}{
The \eqn{\lambda} (scale) parameter controlling the average of the counts.
}
  \item{g}{
The \eqn{\gamma} (shape) parameter controlling variation in the counts.
}
  \item{m}{
The governing parameter of the distribution.
}
}
\details{
This function draws random numbers based on the compound exponential-geomeric series (CEGS) distribution. Data can be fit to this distribution using the \link{cegsB} and \link{cegsML} functions, which respectively use Bayesian and maximum likelihood computations.
The \eqn{\lambda} parameter must be positive, but \eqn{\gamma} may be negative.
Sampling is better when the parameter \eqn{\lambda} is low, not high; variance in the counts is high when \eqn{\gamma} is high in either direction.
The integers returned range from 0 to infinity, so zero counts will probably need to be removed prior to analysing the data.
}
\value{
  A set of randomly generated non-negative integers.
}
\author{
John Alroy
}
\examples{
# draw 100 integers based on a scale parameter of 2 and shape parameter of 4
rcegs(100,l=0.5,g=4)
}
