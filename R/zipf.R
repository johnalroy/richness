zipf<-function(n)	{
	if (length(table(n)) < 3 || max(n) < 3 || is.infinite(max(n[n <= 2^14])))
		return(list('exponent' = NA, 'AICc' = NA, 'fitted.SAD' = NA))
	library(sads)
	z <- coef(fitpower(n))
	p <- dpower(1:2^14,z[1])
	p[p < 1e-15] <- 1e-15
	p <- p / sum(p)
	S <- length(n)
	n2 <- n[n <= 2^14]
	S2 <- length(n2)
	s <- array(dim=2^14,data=0)
	t <- table(n2)
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	ll <- -sum(s[u] * log(p[u]))
	aicc <- 2 * ll + 4 + 4 / (S2 - 2)
	return(list('exponent' = as.numeric(z), 'AICc'= aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
