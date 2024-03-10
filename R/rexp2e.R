rexp2e<-function(n=100,l=1)	{
	sample(0:2^14,n,replace=T,prob=-diff(exp(-l * ((0:(2^14 + 1))^exp(-1)))))
}
