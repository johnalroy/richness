jackknife2<-function(n) {
        N <- sum(n)
        return(length(n) + sum(n == 1) * (2 * N - 3) / N - sum(n == 2) * (N - 2)^2 / (N * (N - 1)))
}
