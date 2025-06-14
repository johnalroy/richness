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

The package fits count vectors to the broken stick (\link{bsd}), compound exponential-geometric series or CEGS (\link{cegsB}, \link{cegsLD}, and \link{cegsML}), discrete Weibull (\link{weibull}), exponential-to-\emph{e} (\link{exp2e}), geometric series (\link{gsd}), half-power exponential (\link{halfpower}), inverse power (\link{ipower}), log series (\link{logd}), negative binomial (\link{nbin}), odds-2 (\link{odds2}), Poisson log normal (\link{pln}), scaled odds (\link{sodds}), zero-sum multinomial (\link{zsm}), and Zipf (\link{zipf}) distributions.
The CEGS, exponential-to-\emph{e}, half-power exponential, inverse power, odds-2, and scaled odds distributions have been overlooked previously.

There are random number generators respectively based on these six distributions and the discrete Weibull (\link{rcegs}, \link{rexp2e}, \link{rhp}, \link{ripower}, \link{rodds2}, \link{rodds}, and \link{rweib}).

The distribution-fitting function \link{saturate} fits a saturated model that assumes a geometric sampling process.

The functions \emph{bsd}, \emph{cegsB}, \emph{cegsLD}, \emph{cegsML}, \emph{exp2e}, \emph{gsd}, \emph{halfpower}, \emph{ipower}, \emph{nbin}, \emph{odds2}, \emph{sodds}, and \emph{weibull} depend only on the utility function \link{sadrad}.
\emph{sadrad} finds a series of counts that are evenly spaced across a species abundance distribution.
The other distribution functions also depend on the packages \emph{poilog} and \emph{sads}.

If a matrix of count data covering a spatial or temporal gradient is available, \link{bat} can be used to carry out before-and-after transect analysis. \emph{bat} can use any one of the many distribution-based richness estimators mentioned previously.

The species richness extrapolation indices are Chao 1 (\link{chao1}), Chao 2 (\link{chao2}), the geometric series index or GSI (\link{gsindex}: Kerr and Alroy 2023), the corrected first-order jackknife (\link{cJ1}: Alroy 2020), the first- and second-order jackknifes (\link{jackknife1} and \link{jackknife2}), and the outer product ratio (\link{opr}).

The blended diversity estimators are Fisher's alpha (\link{fisher}), Shannon's H (\link{shannon}), and the inverse of Simpson's D (\link{simpson}).
\emph{zsm} returns a parameter called theta that is strongly correlated with alpha.

The subsampling methods are rarefaction (\link{rarefy}) shareholder quorum subsampling = coverage-based rarefaction (\link{sqs}: Alroy 2010), and inverted sampling curve extrapolation (\link{isce}).
The latter method provides an asymptotic richness estimate.
It assumes a power law relationship between the inverse of sample size and richness.

The utility functions \link{doubletons} and \link{singletons} can be used to recover key parameters that feature in several estimators.

Most of the distributions are fit by optimising the corrected Akaike information criterion (AICc: Hurvich and Tsai 1993).
The AICc is based on the negative log likelihood of the data.
However, \emph{cegsB} uses a Bayesian calculation and \emph{cegsLD} uses likelihood differencing.

The likelihood is just the product of the probabilities of obtaining the indivual species tallies given the hypothesised species abundance distribution.
Because combinatoric terms are constant across all potential model parameters, this simple calculation yields the same solution as computing multinomial probabilities of the observed count distribution given different parameters.
Packages \emph{poilog} and \emph{sads} use the same likelihood function.

The likelihood expression used in the relevant functions is:

\code{
-sum(s[u] * log(p[u]))
}

where s = the observed SAD (meaning counts of species in each class); p = the hypothesised probabilities of the model; and u = a list of unique counts, i.e., SAD cells that have non-zero counts.

Unlike fitting a function directly to the rank-abundance distribution (RAD), this calculation requires specifying nothing other than the count class probabilities.
So its advantage is being informed strongly by the counts but not beholden to their rank ordering.
The method is therefore not compromised by the fact that the ranks and counts in an RAD are non-independent, which violates the assumptions of fitting methods such as minimising the deviance or sum of squares.

Fitting SADs by likelihood is not robust when the number of singletons is relatively large and the number of occupied classes is small. Therefore, all functions in this package using the AICc fitting approach return no result if the maximum count is one or two (or if there are only one or two classes of counts).

}
\author{
\packageAuthor{richness}

Maintainer: \packageMaintainer{richness}
}
\references{
Alroy, J. 2010. The shifting balance of diversity among major animal groups. \emph{Science} 329:1191-1194.

Alroy, J. 2020. On four measures of taxonomic richness. \emph{Paleobiology} 46:158-175.

Hurvich, C. M. and C. L. Tsai. 1993. A corrected Akaike information criterion for vector autoregressive model selection. \emph{Journal of Time Series Analysis} 14:271-279. 

Kerr, M. R. and J. Alroy. 2023. Body size and abundance are decoupled from species richness in Australian marine bivalves. \emph{Frontiers of Biogeography} 15(4). http://dx.doi.org/10.21425/F5FBG58651
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
# there are 100 species with an average abundance of exp(0) = 1

n <- rpois(100,exp(rnorm(100,sd=2)))
n <- n[n > 0]

# alternatively, create a vector based on the CEGS distribution

#n <- rcegs(1000)
#n <- n[n > 0]

# plot a rank-abundance distribution
# CEGS and the Poisson log normal can look very similar

plot(rev(sort(n)),type='l',log='y',xlab='Rank of count',ylab='Count of individuals')
lines(rev(cegsML(n)$fitted.RAD),col='blue')
lines(rev(logd(n)$fitted.RAD),col='red')
lines(rev(pln(n)$fitted.RAD),col='green4')
legend('topright',c('data','CEGS','LS','PLN'),col=c('black','blue','red','green4'),lwd=1,cex=0.8)

# calculate species richness estimates based on several methods

cegsML(n)$richness
halfpower(n)$richness
ipower(n)$richness
sodds(n)$richness
pln(n)$richness
chao1(n)

# carry out a simulation analysis with 30 samples

# the true number of species R in each sample is log normally distributed
# therefore, the new models may not be accurate (in principle)

R <- ceiling(exp(rnorm(30,mean=log(30))))

# calculate the estimates
# the abundance of each species is drawn with Poisson sampling
# the rates are log-normally distributed
# the rarefaction quota and SQS quorum are low because samples are poor
# some functions may crash if the distributions are extreme

E <- matrix(nrow=30,ncol=22,dimnames=list(1:30,c('R','S','BS','Chao 1','cJ1','GSI','J1','J2','alpha','theta','CEGS','GS','HPE','IP','odds','PLN','Weib','H','D','rare','ISCE','SQS')))

for (i in 1:30)	{
	cat('\rsample',i)
	E[i,1] <- R[i]
	n <- rpois(R[i],exp(rnorm(R[i],sd=2)))
# to simulate the CEGS distribution, use this line instead
#	n <- rcegs(R[i])
	n <- n[n > 0]
	E[i,2] <- length(n)
	if (E[i,2] == 0)
		next
	E[i,3] <- bsd(n)$richness
	E[i,4] <- chao1(n)
	E[i,5] <- cJ1(n)
	E[i,6] <- gsindex(n)
	E[i,7] <- jackknife1(n)
	E[i,8] <- jackknife2(n)
	E[i,9] <- fisher(n)
	E[i,10] <- zsm(n)$theta
	eg <- try(cegsML(n),silent=T)
	if (length(eg) > 1)
		E[i,11] <- eg$richness
	g <- try(gsd(n),silent=T)
	if (length(g) > 1)
		E[i,12] <- g$richness
	E[i,13] <- halfpower(n)$richness
	E[i,14] <- ipower(n)$richness
	E[i,15] <- sodds(n)$richness
	p <- try(pln(n),silent=T)
	if (length(p) > 1)
		E[i,16] <- p$richness
	w <- try(weibull(n),silent=T)
	if (length(w) > 1)
		E[i,17] <- w$richness
	E[i,18] <- shannon(n)
	E[i,19] <- simpson(n)
	E[i,20] <- rarefy(n,25)
	E[i,21] <- isce(n)$asymptotic.richness
	E[i,22] <- sqs(n,0.25)$subsampled.richness / 0.25
}
cat('\n')

# methods that work well yield points that are close to the line of unity

plot(E[,'R'],E[,'CEGS'],xlab='True richness',ylab='Estimate',log='xy',xlim=c(10,1000),ylim=c(10,1000),cex=1.1,pch=21,col='white',bg='blue')
points(E[,'R'],E[,'PLN'],pch=21,col='white',bg='green2')
points(E[,'R'],E[,'SQS'],pch=21,col='white',bg='purple')
points(E[,'R'],E[,'Chao 1'],pch=21,col='white',bg='red')
abline(0,1)
legend(12,1000,legend=c('CEGS','PLN','SQS','Chao 1'),col='white',pt.bg=c('blue','green2','purple','red'),bty='n',pch=21,cex=0.8)

# the ratio of estimated to true richness is close to 1 for a good method
# only methods yielding extrapolative richness estimates are shown here
# the PLN is required to be accurate because the distribution is log normal
# ISCE estimates may be negative if the data are poorly behaved
# therefore, negative estimates are removed

E[E < 0] <- NA

boxplot(E[,c(3:8,11:17,21)] / E[,'R'],cex=0.5,log='y')
abline(0,0)

# ordinate the data using principal coordinates analysis
# estimators that are closer to R (red circle) are more accurate
# the model assumes the PLN, but the PLN is often not accurate

p <- cmdscale(dist(t(log(na.omit(E / E[,'R'])))))
plot(p,cex=0)
text(p,labels=colnames(E),cex=0.7,xpd=NA)
points(p['R',1],p['R',2],col='red',cex=1.5)

# remove measures that reflect evenness strongly and try again
# this does include alpha, theta, rarefaction, and SQS
# also remove ISCE and the Weibull because of their high variance

p <- cmdscale(dist(t(log(na.omit(E[,c(1:8,11:16)] / E[,'R'])))))
plot(p,cex=0)
text(p,labels=colnames(E[,c(1:8,11:16)]),cex=0.7,xpd=NA)
points(p['R',1],p['R',2],col='red',cex=1.5)

# cluster analysis also shows which measures relate to true richness R

plot(hclust(dist(scale(t(na.omit(E))))),ann=F)

}
