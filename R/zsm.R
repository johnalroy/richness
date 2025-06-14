zsm<-function(n)	{
	if (length(table(n)) < 3)
		return(list('J' = NA, 'theta' = NA, 'AICc' = NA, 'fitted.SAD' = NA))
	mx <- max(2^12,2^ceiling(log2(max(n))))
	library(sads)
	z <- coef(fitmzsm(n))
	p <- dmzsm(1:mx,z[1],z[2])
	p[p < 1e-15] <- 1e-15
	p <- p / sum(p)
	S <- length(n)
	s <- array(dim=mx,data=0)
	t <- table(n)
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	ll <- -sum(s[u] * log(p[u]))
	aicc <- 2 * ll + 4 + 12 / (S - 3)
	return(list('J' = as.numeric(z[1]), 'theta' = as.numeric(z[2]), 'AICc'= aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}

