chao2<-function(n,s)	{
	return(length(n) + (s - 1) / s * sum(n == 1) * (sum(n == 1) - 1) / (2 * sum(n == 2) + 1))
}
