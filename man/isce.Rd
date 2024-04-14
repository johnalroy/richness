\name{isce}
\alias{isce}
\title{
Inverted sampling curve extrapolation
}
\description{
Extrapolates a richness estimate by fitting a power law function that relates rarefied richness values to inverted sampling quotas.
}
\usage{
isce(n,base=0.2)
}
\arguments{
  \item{n}{
A vector of integer counts representing the number of individuals per species.
}
  \item{base}{
The lowest sampling proportion used to compute the rarefaction curve.
}
}
\details{
This function fits a power law function to a rarefaction curve, returning a list of parameters that can be used to graph the results.
It first estimates richness at different sampling quotas going down to a \emph{base} proportion of the sample size.
The proportions are equally spaced on an inverse scale.
An example would be 1/5, 1/4, 1/3, 1/2, and 1/1.
So if the sample included 200 individuals, this series would set quotas at 40, 50, 66.7, 100, and 200 individuals.
The quotas are allowed to be fractional values.
The rarefied (= subsampled) richness values are found using the standard combinatorial equation of Hurlbert (1971).

After computing the rarefaction curve, \emph{isce} uses the \emph{stats} function \emph{nls} to fit the power law (1 + \emph{a})/(\emph{p}^\emph{b} + \emph{a}) where \emph{a} = a scaling constant, \emph{p} = is the proportion of individuals drawn in a subsample, and \emph{b} is a power constant.
The data are first log transformed. Asymptotic richness is S (1 + \emph{a})/\emph{a} when sampling is infinite.

This function works well under a variety of circumstances, but can underestimate when the data follow something like the \link{halfpower} distribution. It is generally recommended as a means of illustrating trends in rarefaction curves.
}
\value{
  \item{asymptotic.richness}{The expected number of species based on extrapolating the end of the rarefaction curve.}
  \item{offset}{The parameter governing the scale of the line fitted to the rarefaction curve.}
  \item{power}{The parameter governing the curvature of the line fitted to the rarefaction curve.}
  \item{quotas}{The sampling quotas that were used to compute the rarefaction curve.}
  \item{proportions}{Inverses of the proportions of the total number of individuals that correspond to the quotas.}
  \item{subsampled.richness}{The rarefied species richness values.}
  \item{fitted.curve}{The fitted curve, which extrapolates to the point of infinite sampling.}
}
\references{
Hurlbert, S. H. 1971. The nonconcept of species diversity: a critique and alternative parameters. \emph{Ecology} 52:577-586.
}
\author{
John Alroy
}
\examples{
# Poisson log normal distribution
n <- rpois(100,exp(rnorm(100,sd=2)))
isce(n[n > 0])

# scaled odds distribution
n <- rodds(100,1)
isce(n[n > 0])

# half-power exponential distribution
n <- rhp(100,1)
isce(n[n > 0])

# graph the results
i <- isce(n[n > 0],base=0.1)
plot(i$proportions,i$fitted.curve,type='l',xlim=c(1 / 0.1,0),log='y')
points(i$proportions,i$subsampled.richness,cex=0.5)
}
