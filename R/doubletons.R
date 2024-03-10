doubletons<-function(n,quota)	{
	S <- length(n)
	N <- sum(n)
	if (is.na(N) || N < quota || S == N)
		return(NA)
	all <- lchoose(N , quota)
	s2 <- sum(exp(lchoose(n, 2) + lchoose(N - n, quota - 2) - lchoose(N , quota)))
	return(as.numeric(s2))
}
