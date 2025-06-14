pln<-function(n)	{
	if (length(table(n)) < 3)
		return(list('richness' = NA, 'mu' = NA, 'sigma' = NA, 'AICc' = NA, 'fitted.SAD' = NA))
	S <- length(n)
	mx <- max(2^12,2^ceiling(log2(max(n))))
	pl <- poilog::poilogMLE(n)
	R <- 1000
	q <- seq(1 / (R + 1),1 - 1 / (R + 1),1 / (R + 1))
	le <- qnorm(q,mean=pl$par['mu'],sd=pl$par['sig'])
	e <- exp(le)
	p <- array(dim=mx,data=0)
	dp <- le - e
	p[1] <- mean(exp(dp))
	for (i in 2:mx)	{
		dp <- dp + le - log(i)
		p[i] <- mean(exp(dp))
	}
	p[p < 1e-15] <- 1e-15
	p <- p / sum(p)
	s <- array(dim=mx,data=0)
	t <- table(n[n <= mx])
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	ll <- -sum(s[u] * log(p[u]))
	aicc <- 2 * ll + 4 + 12 / (S - 3)
	return(list('richness' = S / pl$p, 'mu' = as.numeric(pl$par[1]), 'sigma' = as.numeric(pl$par[2]), 'AICc'= aicc, 'fitted.RAD' = sadrad(S,p), 'fitted.SAD' = p[1:2^12]))
}

