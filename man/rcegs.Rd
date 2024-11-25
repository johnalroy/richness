\name{rcegs}
\alias{rcegs}
\title{
Random generation from the CEGS distribution
}
\description{
Draws counts at random from the CEGS distribution.
}
\usage{
rcegs(n=100,l=1,g=0.5)
}
\arguments{
  \item{n}{
The number of counts to be drawn.
}
  \item{l}{
The lambda (scale) parameter controlling the average of the counts.
}
  \item{g}{
The gamma (shape) parameter controlling variation in the counts.
}
  \item{m}{
The governing parameter of the distribution.
}
}
\details{
This function draws random numbers based on the compound exponential-geomeric series (CEGS) distribution. Data can be fit to this distribution using the \link{cegsB} and \link{cegsML} functions, which respectively use Bayesian and maximum likelihood computations.
Both parameters must be positive numbers.
The integers returned range from 0 to infinity, so zero counts will probably need to be removed prior to analysing the data.
}
\value{
  A set of randomly generated non-negative integers.
}
\author{
John Alroy
}
\examples{
# draw 100 integers based on a scale parameter of 2 and shape parameter of 0.25
rcegs(100,l=2,g=0.25)
}
