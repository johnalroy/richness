sqs<-function(n,quorum)	{
	S <- length(n)
	N <- sum(n)
	if (1 - sum(n == 1) / N < quorum)
		return(list('raw.richness' = S, 'individuals' = N, 'subsampled.richness' = NA, 'subsampled.singletons' = NA, 'individuals.drawn' = NA))
	else if (1 - sum(n == 1) / N == quorum)
		return(list('raw.richness' = S, 'individuals' = N, 'subsampled.richness' = S, 'subsampled.singletons' = sum(n == 1), 'individuals.drawn' = N))
	under <- 1
	over <- N
	overS <- 0
	underS <- 0
	overS1 <- 0
	underS1 <- 0
	t <- table(n)
	while (under + 1 < over)	{
		i <- floor(sqrt(under * over)) + 1
		s <- 0
		s1 <- 0
		all <- lchoose(N , i)
		for (j in 1:length(t))	{
			count = as.numeric(names(t[j]))
			s <- s + t[j] - t[j] * exp(lchoose(N - count , i) - all)
			s1 <- s1 + t[j] * exp(log(count) + lchoose(N - count , i - 1) - all)
		}
		if (1 - s1 / i >= quorum)	{
			over <- i
			overS <- s
			overS1 <- s1
		} else	{
			under <- i
			underS <- s
			underS1 <- s1
		}
	}
	return(list('raw.richness' = S, 'individuals' = N, 'subsampled.richness' = as.numeric((underS + overS) / 2), 'subsampled.singletons' = as.numeric((underS1 + overS1) / 2), 'individuals.drawn' = i))
}
