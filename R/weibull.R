weibull<-function(n)	{
	if (length(table(n)) < 3)
		return(list('richness' = NA, 'scale' = NA, 'shape' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	mx <- max(2^12,2^ceiling(log2(max(n))))
	s <- array(dim=mx,data=0)
	t <- table(n)
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	like<-function(a,b)	{
		if (a <= 0 || b <= 0)
			return(1e10)
		p <- (a^u^b - a^(u + 1)^b) / a
		if (is.nan(p[1]) || p[length(p)] < 1e-100 || min(p) == 0)
			return(1e10)
		ll <- -sum(s[u] * log(p))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	cf <- coef(stats4::mle(like,lower=list(a=0,b=0),upper=list(a=1,b=5),start=list(a=0.5,b=0.2)))
	a <- cf[1]
	b <- cf[2]
	if (a == 0 || a ==1)
		return(list('richness' = NA, 'scale' = NA, 'shape' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	aicc <- 2 * like(a,b) + 4 + 12 / (S - 3)
	x <- 1:(mx + 1)
	p <- -diff(a^x^b) / a
	return(list('richness' = as.numeric(S / a), 'scale' = as.numeric(a), 'shape' = as.numeric(b), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
