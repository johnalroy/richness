cegsLD<-function(n)	{
	if (length(table(n)) < 3)
		return(list('richness' = NA, 'scale' = NA, 'shape' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	s <- array(dim=max(n),data=0)
	t <- table(n)
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	px<-function(U,l,g,i)	{
		p <- 1 / ((-log(U))^g / l + 1)
		p[p == 0] <- 1e-4
		dgeom(i,p)
	}
	p0<-function(U,l,g)	{
		1 / ((-log(U))^g / l + 1)
	}
	like<-function(l,g)	{
		p <- array()
		for (i in 1:length(u))
			p[i] <- integrate(px,l=l,g=g,i=u[i],lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
		if (is.nan(p[1]) || p[length(p)] < 1e-100 || min(p) == 0)
			return(1e10)
		p0 <- integrate(p0,l=l,g=g,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
		p <- p / (1 - p0)
		ll <- -sum(s[u] * log(p))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	ml <- matrix(52 / 1:51 - 1,51,51)
	mg <- t(ml)
	ml <- ml / exp(mean(log(n)))
	ll <- matrix(Inf,53,53)
	for (i in 1:51)
		for (j in 1:51)
			ll[i + 1,j + 1] <- like(ml[i,j],mg[i,j])
	ell <- exp(-ll + min(ll) + 5)
	dr <- abs(ell[2:52,2:52] - ell[2:52,1:51]) + abs(ell[2:52,2:52] - ell[2:52,3:53])
	dc <- abs(ell[2:52,2:52] - ell[1:51,2:52]) + abs(ell[2:52,2:52] - ell[3:53,2:52])
	d <- dr * dc
	d <- d / sum(d)
	ml <- 1 / (ml + 1) 
	mg <- 1 / (mg + 1) 
	l <- sum(d * ml)
	g <- sum(d * mg)
	l <- 1 / l - 1
	g <- 1 / g - 1
	aicc <- 2 * like(l,g) + 4 + 12 / (S - 3)
	p <- array()
	mx <- max(2^12,2^ceiling(log2(max(n))))
	for (i in 1:mx)
		p[i] <- integrate(px,l=l,g=g,i=i,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
	p0 <- integrate(p0,l=l,g=g,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
	p <- p / (1 - p0)
	return(list('richness' = as.numeric(S / (1 - p0)), 'scale' = as.numeric(l), 'shape' = as.numeric(g), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}

