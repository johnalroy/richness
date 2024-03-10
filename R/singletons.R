singletons<-function(n,quota)	{
	S <- length(n)
	N <- sum(n)
	if (is.na(N) || N < quota || S == N)
		return(NA)
	s1 <- sum(exp(log(n) + lchoose(N - n, quota - 1) - lchoose(N, quota)))
	return(as.numeric(s1))
}
