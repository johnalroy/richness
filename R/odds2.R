odds2<-function(n)	{
	if (length(n) < 3 || max(n) < 3)
		return(list('richness' = NA, 'mu' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	n2 <- n[n <= 2^14]
	S2 <- length(n2)
	s <- array(dim=2^14,data=0)
	t <- table(n2)
	s[as.numeric(names(t))] <- t
	u <- unique(n2)
	x <- 1:(2^14 + 1)
	like<-function(m)	{
		if (m <= 0)
			return(1e10)
		p <- -diff(1 / (m + x)^2) * (m + 1)^2
		if (is.nan(p[1]) || p[n2[S2]] < 1e-100 || min(p[u]) == 0 || sum(p[u]) > 1 - 1e-120)
			return(1e10)
		ll <- -sum(s[u] * log(p[u]))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	m <- optimise(like,interval=c(0,1e6),maximum=F)$minimum
	if (m == 0 || m == 1e6)
		return(list('richness' = NA, 'mu' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	aicc <- 2 * like(m) + 2 + 4 / (S2 - 2)
	x <- 1:(2^20 + 1)
	p <- -diff(1 / (m + x)^2) * (m + 1)^2
	return(list('richness' = as.numeric(S * (m + 1)^2 / m^2), 'mu' = as.numeric(m), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
