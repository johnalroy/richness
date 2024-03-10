pln<-function(n)	{
	if (length(n) < 3 || max(n) < 3 || is.infinite(max(n[n <= 2^14])))
		return(list('richness' = NA, 'mu' = NA, 'sigma' = NA, 'AICc' = NA, 'fitted.SAD' = NA))
	S <- length(n)
	pl <- poilog::poilogMLE(n)
	R <- 1000
	q <- seq(1 / (R + 1),1 - 1 / (R + 1),1 / (R + 1))
	le <- qnorm(q,mean=pl$par['mu'],sd=pl$par['sig'])
	e <- exp(le)
	p <- array(dim=2^14,data=0)
	dp <- le - e
	p[1] <- mean(exp(dp))
	for (i in 2:2^14)	{
		dp <- dp + le - log(i)
		p[i] <- mean(exp(dp))
	}
	p[p < 1e-15] <- 1e-15
	p <- p / sum(p)
	s <- array(dim=2^14,data=0)
	t <- table(n[n <= 2^14])
	s[as.numeric(names(t))] <- t
	n2 <- n[n <= 2^14]
	S2 <- length(n2)
	u <- unique(n2)
	ll <- -S2 * log(sum(p[u])) - sum(dbinom(s[u],S2,p[u],log=T))
	aicc <- 2 * ll + 4 + 12 / (S2 - 3)
	return(list('richness' = S / pl$p, 'mu' = as.numeric(pl$par[1]), 'sigma' = as.numeric(pl$par[2]), 'AICc'= aicc, 'fitted.RAD' = sadrad(S,p), 'fitted.SAD' = p[1:2^12]))
}
