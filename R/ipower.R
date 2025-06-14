ipower<-function(n)	{
	if (length(table(n)) < 3)
		return(list('richness' = NA, 'beta' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	mx <- max(2^12,2^ceiling(log2(max(n))))
	s <- array(dim=mx,data=0)
	t <- table(n)
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	x <- 2:(mx + 2)
	like<-function(b)	{
		if (b <= 0)
			return(1e10)
		p <- -diff(1 / x^b) * 2^b
		if (is.nan(p[1]) || p[max(n)] < 1e-100 || min(p[u]) == 0 || sum(p[u]) > 1 - 1e-120)
			return(1e10)
		ll <- -sum(s[u] * log(p[u]))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	b <- optimise(like,interval=c(0,10),maximum=F)$minimum
	if (b == 0 || b == 10)
		return(list('richness' = NA, 'beta' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	aicc <- 2 * like(b) + 2 + 4 / (S - 2)
	x <- 2:(mx + 2)
	p <- -diff(1 / x^b) * 2^b
	return(list('richness' = as.numeric(2^b * S), 'beta' = as.numeric(b), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
