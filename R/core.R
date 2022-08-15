bootstrap <- function (n, Bayesian=F, weights=NULL)
{
  if(Bayesian)
  {
    u <- c(0,sort(runif(n-1)),1)
    u <- (u[-1] - u[-length(u)])*n
    if(!is.null(weights)) u <- u*weights*n/sum(u*weights)
  }
  else
  {
    if(is.null(weights)) weights <-rep(1/n,n)
    u <- rmultinom(1,n,weights)
  }
  as.vector(u)
}
