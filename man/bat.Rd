\name{bat}
\alias{bat}
\title{
Before-and-after transect analysis
}
\description{
Calculates species richness and turnover estimates based on a data matrix of counts that represents a spatial or temporal transect.
}
\usage{
bat(x,points,func='sodds')
}
\arguments{
  \item{x}{
A matrix of integer counts giving the number of individuals per species at each point across the transect.
}
  \item{points}{
An optional vector of numbers giving the relative location in space or time of each column in the matrix.
}
  \item{func}{
An optional string specifying the richness estimation function to be used.
}
}
\details{
This function estimates standing richness, boundary richness, species gains, and species losses at points across a spatial or temporal transect.
It can be used with both ecological and palaeontological data.

The method works by estimating the species pool size up to each point across the transect going from left to right ("before" richness or \emph{B}); then again going from right to left ("after" richness or \emph{A}); and then again across the entire data set (total richness or \emph{T}).
Standing richness at each point is then \emph{B} + \emph{A} - \emph{T} based on simple algebra.
Gains are the differences of \emph{B}, losses are the differences of \emph{A}, and boundary richness derives from the same BAT calculation with the points now offset by one step across the transect.

The rows of the data matrix must correspond to species and the columns to sampling points. 
The values must be non-negative integers.

The points must be numbers, but they can be negative. They do not need to be supplied in order, but their ordering must match the column ordering of the matrix.

The computation hinges on comparing different richness estimates, and it only makes sense for them to be computed in the same way.
Therefore, only a single richness estimation function can be supplied.
The richness estimation function must return a "richness" value.
The appropriate functions in this package are \link{bsd}, \link{cegsB}, \link{cegsML}, \link{exp2e}, \link{gsd}, \link{halfpower}, \link{nbin}, \link{odds2}, \link{pln}, \link{sodds}, and \link{weibull}.
\emph{cegsML} is the default because it fits uneven data sets well, but other functions such as \emph{sodds} should be investigated.
The best distribution will yield the lowest AICc values.
Functions not returning richness estimates (such as \link{logd}), returning relative estimates (\link{rarefy} or \link{sqs}), or based on simple indices (such as \link{chao1}) are not designed to work with \emph{bat}.
}
\value{
  \item{total.richness}{The estimated species richness at each point along the transect, including species that only exist at the point.}
  \item{right.edge.richness}{The estimated number of species crossing over between the focal point and the next one to the right.}
  \item{gains}{The estimated number of species gained (first appearing) just before the point.}
  \item{losses}{The estimated number of species lost (last appearing) just after the point.}
}
\author{
John Alroy
}
\examples{
# draw a random matrix of counts from the CEGS distribution
m <- matrix(0,100,10)
n <- rcegs(100,0.25)
for (i in 1:nrow(m))	{
	t <- table(sample(1:10,n[i],replace=T))
	m[i,as.numeric(names(t))] <- t
}
bat(m)

# supply a random sequence of point locations
bat(m,points=runif(10))

# reverse the ordering of the point locations
bat(m,points=10:1)

# use the scaled odds distribution instead of the CEGS distribution
bat(m,func='sodds')
}
