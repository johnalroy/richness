rcegs<-function(n=100,l=1,g=0.5)	{
	rgeom(n,1 / (rexp(n,l) + 1)^(1 / g))
}
