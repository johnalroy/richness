dgeoms<-function()	{
	dg <- matrix(NA,2^14,2^14)
	x <- 1:2^14
	for (i in x)
		dg[i,] <- dgeom(x,1 / (i + 1)) * (i + 1) / i
	dg
}
