\name{richness-package}
\alias{richness-package}
\alias{richness}
\docType{package}
\title{
\packageTitle{richness}
}
\description{
\packageDescription{richness}
}
\details{

This package implements various species abundance distribution, richness, and diversity estimators that work on arrays of species counts.

The package fits count vectors to the broken stick (\link{bsd}), discrete Weibull (\link{weibull}), exponential-to-\emph{e} (\link{exp2e}), geometric series (\link{gsd}), half-power exponential (\link{halfpower}), inverse power (\link{ipower}), log series (\link{logd}), negative binomial (\link{nbin}), odds-2 (\link{odds2}), Poisson log normal (\link{pln}), scaled odds (\link{sodds}), zero-sum multinomial (\link{zsm}), and Zipf (\link{zipf}) distributions.
The exponential-to-\emph{e}, half-power exponential, inverse power, odds-2, and scaled odds distributions have been overlooked previously.
There are random number generators respectively based on these five distributions (\link{rexp2e}, \link{rhp}, \link{ripower}, \link{rodds2}, and \link{rodds}).
The function \link{saturate} fits a saturated model that assumes a geometric sampling process.

The functions \emph{bsd}, \emph{exp2e}, \emph{gsd}, \emph{halfpower}, \emph{ipower}, \emph{nbin}, \emph{odds2}, \emph{sodds}, and \emph{weibull} depend only on the utility function \link{sadrad}.
\emph{sadrad} finds a series of counts that are evenly spaced across a species abundance distribution.
The other distribution functions also depend on the packages \emph{poilog} and \emph{sads}.

If a matrix of count data covering a spatial or temporal gradient is available, \link{bat} can be used to carry out before-and-after transect analysis. \emph{bat} can use any one of the many distribution-based richness estimators mentioned previously.

The species richness extrapolation indices are Chao 1 (\link{chao1}), Chao 2 (\link{chao2}), the corrected first-order jackknife (\link{cJ1}: Alroy 2020), and the first- and second-order jackknifes (\link{jackknife1} and \link{jackknife2}).
The blended diversity estimators are Fisher's alpha (\link{fisher}), Shannon's H (\link{shannon}), and the inverse of Simpson's D (\link{simpson}).
\emph{zsm} returns a parameter called theta that is strongly correlated with alpha.
The subsampling methods are rarefaction (\link{rarefy}) shareholder quorum subsampling = coverage-based rarefaction (\link{sqs}: Alroy 2010), and inverted sampling curve extrapolation (\link{isce}).
The latter method provides an asymptotic richness estimate.
It assumes a power law relationship between the inverse of sample size and richness.

The utility functions \link{doubletons} and \link{singletons} can be used to recover key parameters that feature in several estimators.

Distributions are fit by optimising the corrected Akaike information criterion (AICc: Hurvich and Tsai 1993).
The AICc is based on the negative log likelihood of the data.

Here, the likelihood is just the product of the probabilities of obtaining the indivual species tallies given the hypothesised species abundance distribution.
Because combinatoric terms are constant across all potential model parameters, this simple calculation yields the same solution as computing multinomial probabilities of the observed count distribution given different parameters.
Package \emph{sads} uses the same likelihood function.

The likelihood expression used in the relevant functions is:

\code{
-sum(s[u] * log(p[u]))
}

where s = the observed SAD (meaning counts of species in each class); p = the hypothesised probabilities of the model; and u = a list of unique counts, i.e., SAD cells that have non-zero counts.

Note that only species with 2^14 = 16384 or fewer individuals are considered in each function to make computations tractable.

Unlike fitting a function directly to the rank-abundance distribution (RAD), this calculation requires specifying nothing other than the count class probabilities.
So its advantage is being informed strongly by the counts but not beholden to their rank ordering.
The method is therefore not compromised by the fact that the ranks and counts in an RAD are non-independent, which violates the assumptions of fitting methods such as minimising the deviance or sum of squares.

Fitting SADs by likelihood is not robust when the number of singletons is relatively large and the number of occupied classes is small. Therefore, all functions in this package using the AICc fitting approach return no result if the maximum count is one or two (or if there are only one or two species).

}
\author{
\packageAuthor{richness}

Maintainer: \packageMaintainer{richness}
}
\references{
Alroy, J. 2010. The shifting balance of diversity among major animal groups. \emph{Science} 329:1191-1194.

Alroy, J. 2020. On four measures of taxonomic richness. \emph{Paleobiology} 46:158-175.

Hurvich, C. M. and C. L. Tsai. 1993. A corrected Akaike information criterion for vector autoregressive model selection. \emph{Journal of Time Series Analysis} 14:271-279. 

}
\keyword{ package }
\seealso{
}
\examples{
# please be patient because the examples take some time to run
# some warnings will be printed
# the reason is that two-parameter models can be difficult

# create a single vector of counts (= a sample) using Poisson sampling
# the underlying distribution is log normal and uneven
# there are 100 species with an average abundance of 2

n <- rpois(100,exp(rnorm(100,sd=2)))
n <- n[n > 0]

# alternatively, create a vector based on the scaled odds distribution

n <- rodds(100,1)
n <- n[n > 0]

# plot a rank-abundance distribution
# the Poisson log normal imitates the scaled odds, but not vice versa
# increase the species count to (say) 1000 confirm this

plot(rev(sort(n)),type='l',log='y',xlab='Rank of count',ylab='Count of individuals')
lines(rev(sodds(n)$fitted.RAD),col='blue')
lines(rev(logd(n)$fitted.RAD),col='red')
lines(rev(pln(n)$fitted.RAD),col='green4')
legend('topright',c('data','odds','LS','PLN'),col=c('black','blue','red','green4'),lwd=1,cex=0.8)

# calculate species richness estimates based on several methods

halfpower(n)$richness
ipower(n)$richness
sodds(n)$richness
pln(n)$richness
chao1(n)

# carry out a simulation analysis with 50 samples

# the true number of species R in each sample is log normally distributed
# therefore, the new models are model for the data (in principle)

R <- ceiling(exp(rnorm(50,mean=log(50))))

# calculate the estimates
# the abundance of each species is drawn with Poisson sampling
# the rates are log-normally distributed
# the rarefaction quota and SQS quorum are low because samples are poor
# some functions may crash if the distributions are extreme

E <- matrix(nrow=50,ncol=20,dimnames=list(1:50,c('R','S','BS','Chao 1','cJ1','alpha','theta','GS','J1','J2','HPE','IP','odds','PLN','Weib','H','D','rare','ISCE','SQS')))

for (i in 1:50)	{
	cat('\rsample',i)
	E[i,1] <- R[i]
	n <- rpois(R[i],exp(rnorm(R[i],sd=2)))
# to simulate the scaled odds distribution, use this line instead
#	n <- rodds(R[i])
	n <- n[n > 0]
	E[i,2] <- length(n)
	if (E[i,2] == 0)
		next
	E[i,3] <- bsd(n)$richness
	E[i,4] <- chao1(n)
	E[i,5] <- cJ1(n)
	E[i,6] <- fisher(n)
	E[i,7] <- zsm(n)$theta
	g <- try(gsd(n),silent=T)
	if (length(g) > 1)
		E[i,8] <- g$richness
	E[i,9] <- jackknife1(n)
	E[i,10] <- jackknife2(n)
	E[i,11] <- halfpower(n)$richness
	E[i,12] <- ipower(n)$richness
	E[i,13] <- sodds(n)$richness
	p <- try(pln(n),silent=T)
	if (length(p) > 1)
		E[i,14] <- p$richness
	w <- try(weibull(n),silent=T)
	if (length(w) > 1)
		E[i,15] <- w$richness
	E[i,16] <- shannon(n)
	E[i,17] <- simpson(n)
	E[i,18] <- rarefy(n,25)
	E[i,19] <- isce(n)$asymptotic.richness
	E[i,20] <- sqs(n,0.25)$subsampled.richness / 0.25
}
cat('\n')

# methods that work well yield points that are close to the line of unity

plot(E[,'R'],E[,'odds'],xlab='True richness',ylab='Estimate',log='xy',xlim=c(10,1000),ylim=c(10,1000),cex=1.1,col='blue')
points(E[,'R'],E[,'HPE'],col='gold')
points(E[,'R'],E[,'SQS'],col='purple')
points(E[,'R'],E[,'Chao 1'],col='red')
abline(0,1)
legend(12,1000,legend=c('odds','HPE','SQS','Chao 1'),col=c('blue','gold','purple','red'),bty='n',pch=1,cex=0.7)

# the ratio of estimated to true richness is close to 1 for a good method
# only methods yielding full richness estimates are shown here
# the discrete Weibull is excluded because of its high variance
# the PLN is required to be accurate because the distribution is log normal
# ISCE estimates may be negative if the data are poorly behaved
# therefore, negative estimates are removed

E[E < 0] <- NA

boxplot(E[,c(3:6,8:14,19)] / E[,'R'],log='y')
abline(0,0)

# ordinate the data using principal coordinates analysis
# estimators that are closer to R (red circle) are more accurate

p <- cmdscale(dist(t(log(na.omit(E[,-15] / E[,'R'])))))
plot(p,cex=0)
text(p,labels=colnames(E[,-15]),cex=0.7,xpd=NA)
points(p['R',1],p['R',2],col='red',cex=1.5)

# remove measures that reflect evenness strongly and try again

p <- cmdscale(dist(t(log(na.omit(E[,c(1:14,19)] / E[,'R'])))))
plot(p,cex=0)
text(p,labels=colnames(E[,c(1:14,19)]),cex=0.7,xpd=NA)
points(p['R',1],p['R',2],col='red',cex=1.5)

# cluster analysis also shows which measures relate to true richness

plot(hclust(dist(scale(t(na.omit(E))))),ann=F)

}
