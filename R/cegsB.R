cegsB<-function(n)	{
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
		p <- array()
		for (i in 1:length(u))
			p[i] <- integrate(px,l=l,g=g,i=u[i],lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
		p0 <- integrate(p0,l=l,g=g,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
		p <- p / (1 - p0)
		if (is.na(sum(p)) || is.nan(p[1]) || p[length(p)] < 1e-100 || min(p) == 0)
			return(1e10)
		ll <- -sum(s[u] * log(p))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	q_l <- 1:50 / 51
	pr_l <- matrix(q_l,50,50) / (1 - matrix(q_l,50,50))
	q_g <- 1:25 / 26
	pr_g <- matrix(q_g,50,25,byrow=T) / (1 - matrix(q_g,50,25,byrow=T))
	pr_g <- cbind(-pr_g[,25:1],pr_g)
	ll <- matrix(NA,50,50)
	for (i in 1:50)
		for (j in 1:50)
			ll[i,j] <- like(pr_l[i,j],pr_g[i,j])
	pp <- exp(-(ll - min(ll,na.rm=T)))
	pp <- pp / sum(pp,na.rm=T)
	m_l <- matrix(q_l,50,50)
	m_g <- matrix(q_g,50,25,byrow=T)
	m_g <- cbind(m_g[,25:1],m_g)
	if (sum(pp[,1:25]) > sum(pp[,26:50]))	{
		l <- sum(pp[,1:25] / sum(pp[,1:25]) * m_l[,1:25])
		l <- l / (1 - l)
		g <- sum(pp[,1:25] / sum(pp[,1:25]) * m_g[,1:25])
		g <- -g / (1 - g)
	} else	{
		l <- sum(pp[,26:50] / sum(pp[,26:50]) * m_l[,26:50])
		l <- l / (1 - l)
		g <- sum(pp[,26:50] / sum(pp[,26:50]) * m_g[,26:50])
		g <- g / (1 - g)
	}
	aicc <- 2 * like(l,g) + 4 + 12 / (S2 - 3)
	p <- array()
	for (i in 1:2^14)
		p[i] <- integrate(px,l=l,g=g,i=i,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
	p0 <- integrate(p0,l=l,g=g,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
	p <- p / (1 - p0)
	return(list('richness' = as.numeric(S / (1 - p0)), 'scale' = as.numeric(l), 'shape' = as.numeric(g), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}

