ripower<-function(n=100,b=0.5)	{
	sample(0:2^18,n,replace=T,prob=-diff(1 / (1:(2^18 + 2))^b) * 2^b)
}
