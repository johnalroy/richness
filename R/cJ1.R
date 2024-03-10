cJ1<-function(n)	{
	S <- length(n)
	s1 <- sum(n == 1)
	if (s1 == S)
		return(NA)
	if (s1 == 0)
		return(S)
	l <- log(sum(n) / s1)
	return((S + s1) / (1 - exp(-l) + l * exp(-l)))
}
