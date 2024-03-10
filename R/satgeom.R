satgeom<-function(n,trials=100,dgeoms=NA)	{
	if (length(dgeoms) == 0)
		cat('You must supply a matrix returned by the function dgeoms.\n')
	n <- sort(n)
	S <- length(n)
	n2 <- n[n <= 2^14]
	S2 <- length(n2)
	s <- array(dim=2^14,data=0)
	t <- table(n2)
	s[as.numeric(names(t))] <- t
	u <- unique(n2)
	x <- 1:2^14
	sum_p <- array(dim=2^14,data=0)
	ll <- array(dim=trials)
	for (i in 1:trials)	{
		sub <- n[sample(1:S2,replace=T)]
		s2 <- array(dim=2^14,data=0)
		t <- table(sub)
		s2[as.numeric(names(t))] <- t
		u2 <- unique(sub)
		p <- array(dim=2^14,data=0)
		for (j in u2)
			p <- p + s2[j] * dgeoms[j,]
		p <- p / length(sub)
		sum_p <- sum_p + p
		ll[i] <- -S2 * log(sum(p[u])) - sum(dbinom(s[u],S2,p[u],log=T))
	}
	p <- sum_p / trials
	return(list(mean.LL = mean(ll), LLs = sort(ll), fitted.RAD = sadrad(length(n), p), fitted.SAD = p[1:2^12]))
}
