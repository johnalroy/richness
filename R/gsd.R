gsd<-function(n)	{
	if (length(table(n)) < 3)
		return(list('richness' = NA, 'g' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	mx <- max(2^12,2^ceiling(log2(max(n))))
	s <- array(dim=mx,data=0)
	t <- table(n)
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	x <- 1:mx - 1
	like<-function(g)	{
		g <- g^10
		if (g <= 0 || g >= 1)
			return(1e10)
		p <- (1 - g)^x * g
		if (is.nan(p[1]) || p[max(n)] < 1e-100 || min(p[u]) == 0 || sum(p[u]) > 1 - 1e-120)
			return(1e10)
		ll <- -sum(s[u] * log(p[u]))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	g <- optimise(like,interval=c(0,1)^0.1,maximum=F)$minimum^10
	if (g <= 0 || g >= 1)
		return(list('richness' = NA, 'g' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	aicc <- 2 * like(g^0.1) + 2 + 4 / (S - 2)
	x <- 1:mx - 1
	p <- (1 - g)^x * g
	return(list('richness' = as.numeric(S / (1 - g)), 'p' = as.numeric(g), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
