\name{dgeoms}
\alias{dgeoms}
\title{
Geometric series probability matrix
}
\description{
Computes a matrix of geometric series distribution probabilities based on an array of expected values.
}
\usage{
dgeoms()
}
\details{
This utility function returns a large 2^14 by 2^14 matrix in which each row is a list of probabilities computed by \emph{dgeom} and each column represents a potential parameter of the geometric series distribution.
The parameters correspond to expected integer counts based on the standard geometric series probability mass function. The expectations assume that zero counts will not be observed.
For example, row 3 in the matrix corresponds to an expected value of 3.
This translates into a parameter \emph{p} = 1/(3 + 1), with the density multiplied by 4/3 to make sure the row probabilities sum to 1.
The function is designed to be used in conjunction with \link{satgeom}.
There are no arguments because the size of the matrix is fixed.
\emph{Warning:} this function is very slow, so be patient and be sure to save the output.
}
\author{
John Alroy
}
\examples{
dg <- dgeoms()
satgeom(1:9,100,dg)
}
