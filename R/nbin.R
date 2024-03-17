nbin<-function(n)	{
	if (length(n) < 3 || max(n) < 3)
		return(list('richness' = NA, 'size' = NA, 'probability' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	n2 <- n[n <= 2^14]
	S2 <- length(n2)
	s <- array(dim=2^14,data=0)
	t <- table(n2)
	s[as.numeric(names(t))] <- t
	u <- unique(n2)
	k <- 1:2^14
	like<-function(r,p)	{
		if (r < 1e-4 || r >= 100 || p < 1e-4 || p >= 1)
			return(1e6)
		pr <- exp(lfactorial(k + r) - lfactorial(k + 1) - lfactorial(r - 1) + (k + 1) * log(1 - p) + r * log(p))
		pr <- pr / (1 - p^r)
		if (is.infinite(pr[1]) || sum(pr) == 0)
			return(1e6)
		ll <- -sum(s[u] * log(pr[u]))
		if (is.infinite(ll) || is.nan(ll))
			ll <- 1e6
		ll
	}
	c <- coef(stats4::mle(like,lower=list(r=1e-4,p=1e-4),upper=list(r=100,p=1),start=list(r=1,p=0.5)))
	r <- c[1]
	p <- c[2]
	if (r == 1e-4 || r == 100 || p == 1e-4 || p == 1)
		return(list('richness' = NA, 'size' = NA, 'probability' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	pr <- exp(lfactorial(k + r) - lfactorial(k + 1) - lfactorial(r - 1) + (k + 1) * log(1 - p) + r * log(p))
	pr <- pr / (1 - p^r)
	aicc <- 2 * like(r,p) + 4 + 12 / (S2 - 3)
	return(list('richness' = as.numeric(S / (1 - p^r)), 'size' = as.numeric(r), 'probability' = as.numeric(p), 'AICc' = as.numeric(aicc), 'fitted.RAD' = sadrad(S,pr), 'fitted.SAD' = as.numeric(pr[1:2^12])))
}
