cegsML<-function(n)	{
	if (length(table(n)) < 3 || max(n) < 3)
		return(list('richness' = NA, 'scale' = NA, 'shape' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	n2 <- n[n <= 2^14]
	S2 <- length(n2)
	s <- array(dim=max(n2),data=0)
	t <- table(n2)
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
		if (l <= 0 || is.na(g))
			return(1e10)
		p <- array()
		for (i in 1:length(u))
			p[i] <- integrate(px,l=l,g=g,i=u[i],lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
		p0 <- integrate(p0,l=l,g=g,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
		p <- p / (1 - p0)
		if (is.nan(p[1]) || p[length(p)] < 1e-100 || min(p) == 0)
			return(1e10)
		ll <- -sum(s[u] * log(p))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	cf <- try(coef(stats4::mle(like,lower=list(l=0,g=-1),upper=list(l=1e8,g=10),start=list(l=1,g=2))),silent=T)
	cf2 <- try(coef(stats4::mle(like,lower=list(l=0,g=-10),upper=list(l=1e8,g=1),start=list(l=1,g=-2))),silent=T)
	if (length(cf) == 1 && length(cf2) == 1)
		return(list('richness' = NA, 'scale' = NA, 'shape' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	if (length(cf) == 1 || like(cf2[1],cf2[2]) < like(cf[1],cf[2]))
		cf <- cf2
	l <- cf[1]
	g <- cf[2]
	if (l == 0 || l == 1e8 || g == -10 || g == 10)
		return(list('richness' = NA, 'scale' = NA, 'shape' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	aicc <- 2 * like(l,g) + 4 + 12 / (S2 - 3)
	p <- array()
	for (i in 1:2^14)
		p[i] <- integrate(px,l=l,g=g,i=i,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
	p0 <- integrate(p0,l=l,g=g,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
	p <- p / (1 - p0)
	return(list('richness' = as.numeric(S / (1 - p0)), 'scale' = as.numeric(l), 'shape' = as.numeric(g), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
