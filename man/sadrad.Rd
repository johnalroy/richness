\name{sadrad}
\alias{sadrad}
\title{
Predicted rank abundance distributions
}
\description{
Computes a rank abundance distribution from a species abundance distribution.
}
\usage{
sadrad(S,p)
}
\arguments{
  \item{s}{
An integer representing the number of individuals in the sample.
}
  \item{p}{
A vector including the predicted proportions of species belonging to each abundance class, starting with the singletons and doubletons.
}
}
\details{
This utility function computes predicted integer counts for a rank abundance distribution (RAD) based on a specified species abundance probability distribution (SAD). The algorithm is used throughout the \emph{richness} package.
}
\author{
John Alroy
}
\examples{
# compute an artificial SAD with 1024 count categories
p <- sort(exp(rnorm(1024,sd=3)),decreasing=T)
p <- p / sum(p)

# predict the RAD for 100 species based on the random distribution
sadrad(100,p)
}
