rodds<-function(n=100,m=1)	{
	x <- 0:2^18
	sample(x,n,replace=T,prob=m / ((m + x) * (m + x + 1)))
}
