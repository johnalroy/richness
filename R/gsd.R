gsd<-function(n)	{
	if (length(n) < 3 || max(n) < 3)
		return(list('richness' = NA, 'k' = NA, 'lambda' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	n2 <- n[n <= 2^14]
	S2 <- length(n2)
	s <- array(dim=2^14,data=0)
	t <- table(n2)
	s[as.numeric(names(t))] <- t
	u <- unique(n2)
	x <- 1:2^14
	x1 <- x + 1
	like<-function(l)	{
		l <- exp(-l)
		if (l <= -100)
			return(1e10)
		p <- exp(x * log(1 / l) - x1 * log(1 / l + 1))
		p <- p / sum(p)
		if (is.nan(p[1]) || p[n2[S2]] < 1e-100 || min(p[u]) == 0 || sum(p[u]) > 1 - 1e-120)
			return(1e10)
		ll <- -S2 * log(sum(p[u])) - sum(dbinom(s[u],S2,p[u],log=T))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	l <- coef(stats4::mle(like,lower=list(l=-100),upper=list(l=1e3),start=list(l=log(mean(n)))))
	if (l == -100 || l == 1e3)
		return(list('richness' = NA, 'k' = NA, 'lambda' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	aicc <- 2 * like(l) + 2 + 4 / (S2 - 2)
	x <- 1:2^20
	x1 <- x + 1
	l <- exp(-l)
	p <- exp(x * log(1 / l) - x1 * log(1 / l + 1)) * (l + 1)
	return(list('richness' = as.numeric(S * (l + 1)), 'k' = as.numeric(l / (l + 1)), 'lambda' = as.numeric(l), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
