bsd<-function(n)	{
	if (length(n) < 3 || max(n) < 3)
		return(list('richness' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	R <- length(n)
	n2 <- n[n <= 2^14]
	R2 <- length(n2)
	N <- sum(n2)
	s <- array(dim=2^14,data=0)
	t <- table(n2)
	s[as.numeric(names(t))] <- t
	u <- unique(n2)
	like<-function(S)	{
		p <- (1 - 1:N / N)^(S - 2) * (S - 1) / N
		p <- p[1:N] / sum(p[1:N])
		if (length(p) < 2^14)
			p[(N + 1):2^14] <- 0
		else
			p <- p[1:2^14]
		if (min(p[u]) == 0 || sum(p[u]) > 1 - 1e-120)
			return(1e10)
		ll <- -sum(s[u] * log(p[u]))
		if (is.infinite(ll))
			return(1e10)
		ll
	}
	S <- coef(stats4::mle(like,lower=list(S=length(n)),upper=list(S=100 * length(n)),start=list(S=2 * length(n))))
	if (S == 100 * length(n))	{
		warning('richness estimate is out of range')
		return(list('richness' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	}
	p <- (1 - 1:N / N)^(S - 2) * (S - 1) / N
	if (length(p) < 2^14)
		p[(N + 1):2^14] <- 0
	else
		p <- p[1:2^14]
	aicc <- 2 * like(S) + 2 + 4 / (length(n2) - 2)
	return(list('richness' = as.numeric(S), 'AICc' = as.numeric(aicc), 'fitted.RAD' = sadrad(R,p), 'fitted.SAD' = p[1:2^12]))
}
