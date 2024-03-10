weibull<-function(n)	{
	if (length(n) < 3 || max(n) < 3)
		return(list('richness' = NA, 'scale' = NA, 'shape' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	n2 <- n[n <= 2^14]
	S2 <- length(n2)
	s <- array(dim=2^14,data=0)
	t <- table(n2)
	s[as.numeric(names(t))] <- t
	u <- unique(n2)
	x <- 1:2^14
	x1 <- 1:2^14 + 1
	like<-function(a,b)	{
		if (a < 1e-6 || b < 1e-6)
			return(1e10)
		p <- (exp(-(x / a)^b) - exp(-(x1 / a)^b)) / exp(-1 / a^b)
		if (is.nan(p[1]) || sum(p) == 0 || min(p[u]) == 0 || sum(p[u]) > 1 - 1e-120)
			return(1e10)
		ll <- -S2 * log(sum(p[u])) - sum(dbinom(s[u],S2,p[u],log=T))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	c <- coef(stats4::mle(like,lower=list(a=1e-6,b=1e-6),upper=list(a=1e2,b=5),start=list(a=1,b=0.2)))
	a <- c[1]
	b <- c[2]
	if (a < 1e-5 || a > 1e3 - 1 || b < 1e-5 || b > 9.99)
		return(list('richness' = NA, 'scale' = NA, 'shape' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	p <- (exp(-(x / a)^b) - exp(-(x1 / a)^b)) / exp(-1 / a^b)
	aicc <- 2 * like(a,b) + 4 + 12 / (S2- 3)
	return(list('richness' = as.numeric(S / exp(-1 / a^b)), 'scale' = as.numeric(a), 'shape' = as.numeric(b), 'AICc' = as.numeric(aicc), 'fitted.RAD' = sadrad(S,p), 'fitted.SAD' = as.numeric(p[1:2^12])))
}
