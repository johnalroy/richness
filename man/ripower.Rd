\name{ripower}
\alias{ripower}
\title{
Random generation from the inverse power distribution
}
\description{
Draws counts at random from the inverse power distribution.
}
\usage{
ripower(n=100,b=0.5)
}
\arguments{
  \item{n}{
The number of counts to be drawn.
}
  \item{l}{
The governing parameter of the distribution (\emph{beta}).
}
}
\details{
This function draws random numbers based on the inverse power distribution. Data can be fit to this distribution using the \link{ipower} function.
The governing parameter l (\emph{beta}) is a measure of sampling intensity. A value of around 0.5 is commonly encountered in real ecological data.
The integers returned range from 0 to 2^14 (= 16384), so zero counts may need to be removed prior to analysing the data.
}
\value{
  A set of randomly generated non-negative integers.
}
\author{
John Alroy
}
\examples{
# draw 100 integers based on a beta parameter of 0.5
ripower(100,0.5)
}
