\name{rhp}
\alias{rhp}
\title{
Random generation from the half-power exponential distribution
}
\description{
Draws counts at random from the half-power exponential distribution.
}
\usage{
rhp(n=100,l=1)
}
\arguments{
  \item{n}{
The number of counts to be drawn.
}
  \item{l}{
The governing parameter of the distribution (\emph{lambda}).
}
}
\details{
This function draws random numbers based on the half-power exponential distribution. Data can be fit to this distribution using the \link{halfpower} function.
The governing parameter l (\emph{lambda}) is a measure of sampling intensity. A value of around 1 is commonly encountered in real ecological data.
The integers returned range from 0 to 2^14 (= 16384), so zero counts will probably need to be removed prior to analysing the data.
}
\value{
  A set of randomly generated non-negative integers.
}
\author{
John Alroy
}
\examples{
# draw 100 integers based on a lambda parameter of 1
rhp(100,1)
}
