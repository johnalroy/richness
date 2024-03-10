chao1<-function(n)	{
	return(length(n) + sum(n == 1) * (sum(n == 1) - 1) / (2 *  (sum(n == 2) + 1)))
}
