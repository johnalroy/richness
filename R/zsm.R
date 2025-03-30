zsm<-function(n)	{
	if (length(table(n)) < 3 || max(n) < 3 || is.infinite(max(n[n <= 2^14])))
		return(list('J' = NA, 'theta' = NA, 'AICc' = NA, 'fitted.SAD' = NA))
	library(sads)
	z <- coef(fitmzsm(n))
	p <- dmzsm(1:2^14,z[1],z[2])
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
	aicc <- 2 * ll + 4 + 12 / (S2 - 3)
	return(list('J' = as.numeric(z[1]), 'theta' = as.numeric(z[2]), 'AICc'= aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
