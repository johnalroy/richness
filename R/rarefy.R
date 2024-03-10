rarefy<-function(n,quota)	{
	N <- sum(n)
	if (is.na(N) || N < quota || length(n) == N)
		return(NA)
	r <- sum(1 - exp(lchoose(N - n, quota) - lchoose(N, quota)))
	return(as.numeric(r))
}
