exp2e<-function(n)	{
	if (length(table(n)) < 3)
		return(list('richness' = NA, 'lambda' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	mx <- max(2^12,2^ceiling(log2(max(n))))
	s <- array(dim=mx,data=0)
	t <- table(n)
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	e <- exp(-1)
	x <- (1:(mx + 1))^e
	like<-function(l)	{
		if (l <= -1)
			return(1e10)
		p <- -diff(exp(-l * x)) / exp(-l)
		if (is.nan(p[1]) || p[max(n)] < 1e-100 || min(p[u]) == 0 || sum(p[u]) > 1 - 1e-120)
			return(1e10)
		ll <- -sum(s[u] * log(p[u]))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	l <- optimise(like,interval=c(-1,10),maximum=F)$minimum
	if (l == -1 || l == 100)
		return(list('richness' = NA, 'lambda' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	aicc <- 2 * like(l) + 2 + 4 / (S - 2)
	x <- (1:(mx + 1))^e
	p <- -diff(exp(-l * x)) / exp(-l)
	return(list('richness' = as.numeric(S / e^l), 'lambda' = as.numeric(l), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
