logd<-function(n)	{
	if (length(table(n)) < 3)
		return(list('alpha' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	a <- fisher(n)
	S <- length(n)
	N <- sum(n)
	mx <- max(2^12,2^ceiling(log2(max(n))))
	x <- N / (N + a)
	p <- -1 / log(1 - x) * x^(1:mx) / 1:mx
	p[p < 1e-15] <- 1e-15
	s <- array(dim=mx,data=0)
	t <- table(n)
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	ll <- -sum(s[u] * log(p[u]))
	aicc <- 2 * ll + 2 + 4 / (S - 2)
	return(list('alpha' = as.numeric(a), 'AICc'= aicc, 'fitted.RAD' = sadrad(S,p), 'fitted.SAD' = p[1:2^12]))
}

