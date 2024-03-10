isce<-function(n,base=0.2)	{
	library(minpack.lm)
	library(stats4)
	N <- sum(n)
	S <- length(n)
	if (S < 4)
		return(list(asymptotic.richness = NA, offset = NA, power = NA, quotas = NA, proportions = NA, subsampled.richness = NA, fitted.curve = NA))
	lc<-function(a,b)	{
		lfactorial(a) - lfactorial(a - b) - lfactorial(b)
	}
	rarefy<-function(q)	{
		sum(1 - exp(lc(N - n, q) - lc(N, q)))
	}
	p <- 0:100 / (base * 100)
	q <- N / p
	s <- array()
	for (i in which(p == 1):length(p))
		s[i] <- rarefy(q[i])
	q[p < 1] <- NA
	q1 <- q[p >= 1]
	s1 <- s[p >= 1]
	p1 <- p[p >= 1]
	cf <- coef(nlsLM(log(s1) ~ log(S * (1 + a) / (p1^b + a)),start=list(a=1,b=1)))
	a <- cf[1]
	b <- cf[2]
	list(asymptotic.richness = as.numeric(S * (1 + a) / a), offset = as.numeric(a), power = as.numeric(b), quotas = q, proportions = p, subsampled.richness = s, fitted.curve = S * (1 + a) / (p^b + a))
}
