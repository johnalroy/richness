bat<-function(x,points,func='sodds')	{
	if (! is.matrix(x) || ! is.numeric(x))	{
		warning('The first supplied argument must be a matrix of non-zero integers')
		return(NA)
	}
	if (missing(points))
		points <- 1:ncol(x)
	x <- x[rowSums(x) > 0,]
	points <- points[colSums(x) > 0]
	x <- x[,colSums(x) > 0]
	u <- unique(points)
	l <- length(u)
	b <- array()
	a <- array()
	for (i in 1:l)	{
		n <- x[,points <= u[i]]
		if (is.matrix(n))
			n <- rowSums(n)
		b[i] <- do.call(func,list(n[n > 0]))$richness
		n <- x[,points >= u[i]]
		if (is.matrix(n))
			n <- rowSums(n)
		a[i] <- do.call(func,list(n[n > 0]))$richness
	}

	list(total.richness = b + a - do.call(func,list(rowSums(x)))$richness,
	right.edge.richness = c(b[1:(l - 1)] + a[2:l] - do.call(func,list(rowSums(x)))$richness,NA),
	gains = c(NA,diff(b)),
	losses = c(diff(-a),NA))
}
