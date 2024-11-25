\name{rweib}
\alias{rweib}
\title{
Random generation from the discrete Weibull distribution
}
\description{
Draws counts at random from the discrete Weibull distribution.
}
\usage{
rweib(n=100,a=0.5,b=0.5)
}
\arguments{
  \item{n}{
The number of counts to be drawn.
}
  \item{a}{
The scale parameter controlling the average of the counts.
}
  \item{b}{
The shape parameter controlling variation in the counts.
}
  \item{m}{
The governing parameter of the distribution.
}
}
\details{
This function draws random numbers based on the discrete Weibull distribution. Data can be fit to this distribution using the \link{weibull} function.
The scale parameter must be greater than 0 and less than 1.
The shape parameter must be a positive number.
The integers returned range from 0 to 2^18 (= 262144), so zero counts will probably need to be removed prior to analysing the data.
}
\value{
  A set of randomly generated non-negative integers.
}
\author{
John Alroy
}
\examples{
# draw 100 integers based on a scale parameter of 0.9 and shape parameter of 1
rweib(100,a=0.9,b=1)
}
