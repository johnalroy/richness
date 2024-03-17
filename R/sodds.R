sodds<-function(n)	{
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
	x <- 1:2^14
	like<-function(m)	{
		if (m <= -1)
			return(1e10)
		p <- (m + 1) / ((m + x) * (m + x + 1))
		p <- p / sum(p)
		if (is.nan(p[1]) || p[n2[S2]] < 1e-100 || min(p[u]) == 0 || sum(p[u]) > 1 - 1e-120)
			return(1e10)
		ll <- -sum(s[u] * log(p[u]))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	m <- coef(stats4::mle(like,lower=list(m=-1),upper=list(m=1e4),start=list(m=1)))
	if (m == -1 || m == 1e4)
		return(list('richness' = NA, 'mu' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	aicc <- 2 * like(m) + 2 + 4 / (S2 - 2)
	x <- 1:2^20
	p <- (m + 1) / ((m + x) * (m + x + 1))
	p <- p / sum(p)
	r <- NA
	if (m > 0 && exp(mean(log(n))) > 2.2 && s[1] / length(n) < 0.5)
		r <- S * (m + 1) / m
	return(list('richness' = as.numeric(r), 'mu' = as.numeric(m), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
