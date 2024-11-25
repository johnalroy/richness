\name{rexp2e}
\alias{rexp2e}
\title{
Random generation from the exponential-to-e distribution
}
\description{
Draws counts at random from the exponential-to-\emph{e} (exp2e) distribution.
}
\usage{
rexp2e(n=100,l=1)
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
This function draws random numbers based on the exponential-to-\emph{e} distribution. Data can be fit to this distribution using the \link{exp2e} function.
The governing parameter l (\emph{lambda}) is a measure of sampling intensity. A value of around 1 is commonly encountered in real ecological data.
The integers returned range from 0 to 2^14 (= 16384), so zero counts may need to be removed prior to analysing the data.
}
\value{
  A set of randomly generated non-negative integers.
}
\author{
John Alroy
}
\examples{
# draw 100 integers based on a lambda parameter of 1
rexp2e(100,1)
}
