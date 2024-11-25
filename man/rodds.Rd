\name{rodds}
\alias{rodds}
\title{
Random generation from the scaled odds distribution
}
\description{
Draws counts at random from the scaled odds distribution.
}
\usage{
rodds(n=100,m=1)
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
This function draws random numbers based on the scaled odds distribution. Data can be fit to this distribution using the \link{sodds} function.
The governing parameter m (\emph{mu}) is literally the median expected number of individuals across all species in the community (not only the known species). A value of around 4 is commonly encountered in real ecological data.
The integers returned range from 0 to 2^14 (= 16384), so zero counts will probably need to be removed prior to analysing the data.
}
\value{
  A set of randomly generated non-negative integers.
}
\author{
John Alroy
}
\examples{
# draw 100 integers based on a mu parameter of 1
rodds(100,1)
}
