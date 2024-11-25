rweib<-function(n=100,a=0.5,b=0.5)	{
	sample(0:2^18,n,replace=T,prob=-diff(a^(0:(2^18 + 1))^b))
}
