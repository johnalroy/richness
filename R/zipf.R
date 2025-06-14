zipf<-function(n)	{
	if (length(table(n)) < 3)
		return(list('exponent' = NA, 'AICc' = NA, 'fitted.SAD' = NA))
	library(sads)
	z <- coef(fitpower(n))
	mx <- max(2^12,2^ceiling(log2(max(n))))
	p <- dpower(1:mx,z[1])
	p[p < 1e-15] <- 1e-15
	p <- p / sum(p)
	S <- length(n)
	s <- array(dim=mx,data=0)
	t <- table(n)
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	ll <- -sum(s[u] * log(p[u]))
	aicc <- 2 * ll + 4 + 4 / (S - 2)
	return(list('exponent' = as.numeric(z), 'AICc'= aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
