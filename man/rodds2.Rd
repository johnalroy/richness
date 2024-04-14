\name{rodds2}
\alias{rodds2}
\title{
Random generation from the odds-2 distribution
}
\description{
Draws counts at random from the odds-2 series distribution.
}
\usage{
rodds2(n=100,m=3)
}
\arguments{
  \item{n}{
The number of counts to be drawn.
}
  \item{m}{
The governing parameter of the distribution.
}
}
\details{
This function draws random numbers based on the odds-2 series distribution. Data can be fit to this distribution using the \link{odds2} function.
The governing parameter is m (\emph{mu}).
The median expected number of individuals across all species in the community (not only the known species) is m/0.5^0.5 - m, so if m = 10 the expected median is 4.14.
A value above 10 is commonly encountered in real ecological data, so the default value is low.
The integers returned range from 0 to 2^14 (= 16384), so zero counts may need to be removed prior to analysing the data.
}
\value{
  A set of randomly generated non-negative integers.
}
\author{
John Alroy
}
\examples{
# draw 100 integers based on a mu parameter of 3
rodds2(100,3)
}
