jackknife1<-function(n) {
        N <- sum(n)
        return(length(n) + sum(n == 1) * (N - 1) / N)
}
