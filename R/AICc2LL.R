AICc2LL<-function(n,AICc,k=1)	{
	(AICc - 2 * k - (2 * k^2 + 2 * k) / (sum(n <= 2^14) - k - 1)) / 2
}
