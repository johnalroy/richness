opr<-function(n)	{
	if (length(n) == 0 || max(n) < 2)
		return(NA)
	length(n) / mean((n %o% (n - 1))^0.5 / (n %o% n)^0.5)
}
