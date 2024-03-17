logd<-function(n)	{
	if (length(n) < 3 || max(n) < 3 || is.infinite(max(n[n <= 2^14])))
		return(list('alpha' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	a <- fisher(n)
	S <- length(n)
	N <- sum(n)
	x <- N / (N + a)
	p <- -1 / log(1 - x) * x^(1:2^16) / 1:2^16
	p[p < 1e-15] <- 1e-15
	n2 <- n[n <= 2^14]
	S2 <- length(n2)
	s <- array(dim=2^14,data=0)
	t <- table(n2)
	s[as.numeric(names(t))] <- t
	u <- unique(n2)
	ll <- -sum(s[u] * log(p[u]))
	aicc <- 2 * ll + 2 + 4 / (S2 - 2)
	return(list('alpha' = as.numeric(a), 'AICc'= aicc, 'fitted.RAD' = sadrad(S,p), 'fitted.SAD' = p[1:2^12]))
}
