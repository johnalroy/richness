rodds2<-function(n,m=3)	{
	sample(0:2^18,n,replace=T,prob=-diff(1 / (m + 0:(2^18 + 1))^2) * (m + 1)^2)
}
