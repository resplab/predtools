#' @title EVPI (Expected Value of Perfect Information) for validation
#' Takes a vector of mean and a 2X2 covariance matrix
#' @param Y Binary response variable
#' @param pi Mean of the second distribution
#' @param method EVPI calculation method
#' @param n_sim Number of Monte Carlo simulations (for bootstrap-based methods)
#' @param zs vector of risk thresholds at which EVPI is to be calculated
#' @param weights (optional) observation weights 
#' @return Returns a data frame containing thresholds, EVPIs, and some auxilary output.
#' @export
evpi_val <- function(Y, pi, method=c("bootstrap","bayesian_bootstrap","asymptotic"), n_sim=1000, zs=(0:99)/100, weights=NULL)
{
  n <- length(Y)
  
  if(method=="asymptotic")
  {
    if(is.null(weights)) weights <- rep(1,n)
    
    ENB_perfect <- ENB_current <- rep(0, length(zs))
    for(j in 1:length(zs))
    {
      NB_model <- sum(weights*(pi>zs[j])*(Y-(1-Y)*zs[j]/(1-zs[j])))/n
      NB_all <- sum(weights*(Y-(1-Y)*zs[j]/(1-zs[j])))/n
      parms <- calc_NB_moments(Y,pi,zs[j], weights)
      if(is.na(parms[5])){
        ENB_perfect[j] <- ENB_current[j] <- max(0,NB_model,NB_all)
      } else{
        if(parms[5]>0.999999) parms[5]<- 0.999999
        if(parms[5]< -0.999999) parms[5]<- -0.999999
        
        
        tryCatch(
          {ENB_perfect[j] <- mu_max_trunc_bvn(parms[1],parms[2],parms[3],parms[4],parms[5])}, error=function(cond) {
            return(NULL)
          })
        ENB_current[j] <- max(0,NB_model,NB_all)
      }
    }
    return(data.frame(z=zs, ENB_perfect=ENB_perfect, ENB_current=ENB_current, EVPIv=ENB_perfect-ENB_current))
  }
  
  NB_model <- NB_all <- matrix(0, n_sim, ncol=length(zs))
  
  if(method=="bootstrap" || method=="bayesian_bootstrap")
  {
    Bayesian_bootstrap <- method=="bayesian_bootstrap"
    for(i in 1:n_sim)
    {
      w_x <- bootstrap(n, Bayesian_bootstrap, weights = weights)
      for(j in 1:length(zs))
      {
        NB_model[i,j] <- sum(w_x*(pi>zs[j])*(Y-(1-Y)*zs[j]/(1-zs[j])))/n
        NB_all[i,j] <- sum(w_x*(Y-(1-Y)*zs[j]/(1-zs[j])))/n
      }
    }
  }
  else
  {
    stop("Method ",method," is not recognized.")
  }
  
  ENB_model <- ENB_all <- ENB_perfect <- ENB_current <- EVPIv <- p_useful <- rep(NA,length(zs))
  for(i in 1:length(zs))
  {
    ENB_model[i] <- mean(NB_model[,i])
    ENB_all[i] <- mean(NB_all[,i])
    
    ENB_perfect[i] <- mean(pmax(NB_model[,i],NB_all[,i],0))
    ENB_current[i] <- max(ENB_model[i],ENB_all[i],0)
    EVPIv[i] <- ENB_perfect[i] - ENB_current[i]
    p_useful[i] <- mean((pmax(NB_model[,i],NB_all[,i],0)-NB_model[,i])==0)
  }
  
  data.frame(z=zs, ENB_model=ENB_model, ENB_all=ENB_all, ENB_current=ENB_current, ENB_perfect=ENB_perfect, EVPIv=EVPIv, p_useful=p_useful)
}





#' @title Calculates the expected value of the maximum of two random variables with zero-truncated bivariate normal distribution
#' Takes a vector of mean and a 2X2 covariance matrix
#' @param mu1 Mean of the first distribution
#' @param mu2 Mean of the second distribution
#' @param sigma1 SD of the first distribution
#' @param sigma2 SD of the second distribution
#' @param rho Correlation coefficient of the two random variables
#' @param precision Numerical precision value
#' @return A scalar value for the expected value
#' @export
mu_max_trunc_bvn <-function(mu1, mu2, sigma1, sigma2, rho, precision=.Machine$double.eps) {
  
  # indicator function:
  # returns 1 if the condition is true
  #         0 otherwise
  indicator <- function(condition){
    as.numeric(condition)
  }
  
  # To deal with numerical precision,
  # we test whether abs(Q) < precision
  # instead of testing for Q=0
  # Machine precision is 2.220446e-16
  
  # the first term of the integral:
  # uv evaluated at infinity - uv evaluated at 0
  uv <-  function(mu1, mu2, sigma1, sigma2, rho) {
    # two commonly used quantities
    tmp1 <- sigma1 - rho * sigma2
    tmp2 <- (rho * sigma2 * mu1 - sigma1 * mu2 ) / (sigma1 * sigma2 * sqrt(1 - rho ^2))
    
    # compute [uv|_0^\infinity
    mu1 * (indicator(tmp1 > precision) +  indicator(abs(tmp1) < precision) * pnorm(tmp2)) -
      pnorm(tmp2) * (-sigma1 * dnorm(-mu1 / sigma1) + mu1 * pnorm(-mu1 / sigma1))
  }
  
  # the second term of the integral:
  # integral of vdu from 0 to infinity
  vdu <- function(mu1, mu2, sigma1, sigma2, rho) {
    # quantities used many times
    tmp1 <- (sigma1 - rho * sigma2)
    delta <- sign(tmp1)
    alpha <- (sigma1 * mu2 - rho * sigma2 * mu1)
    beta <- sigma1 * sigma2 * sqrt(1 - rho ^ 2)
    
    # the first case:
    # sigma1 - rho * sigma2 = 0
    if(abs(tmp1)<precision){
      return(0)
    } 
    # the second case:
    # sigma1 - rho * sigma2 \not = 0
    else{
      
      # a term used multiple times as the denominator
      denom <- sqrt(sigma1^2+sigma2^2-2*rho*sigma1*sigma2)
      
      # correlation for the standard bivariate normal term
      sBVN_rho <- -sigma2*sqrt(1-rho^2)/denom
      
      # first term T1
      T1 <- delta*mu1*(pnorm(delta*(mu2-mu1)/denom) - as.numeric(pmvnorm(
        lower = c(-Inf,-Inf),
        upper = c(delta*(mu2-mu1)/denom,-delta*alpha / beta),
        mean =
          c(0, 0),
        sigma = matrix(c(1, sBVN_rho, sBVN_rho, 1), 2)
      )[[1]]) )
      
      # quantities for T2
      a <- (mu2-mu1)/tmp1
      b <- (beta/abs(tmp1))/sigma1
      t <- sqrt(1 + b ^ 2)
      
      # second term T2
      T2 <- -delta*sigma1 / t * dnorm(a / t) *
        (1 - pnorm(-delta * t * alpha / beta + a * b / t))
      
      return(T1+T2)
    }
  }
  
  # compute E(max(Y,0))
  uv(mu1, mu2, sigma1, sigma2, rho) + uv(mu2, mu1, sigma2, sigma1, rho) -
    vdu(mu1, mu2, sigma1, sigma2, rho) - vdu(mu2, mu1, sigma2, sigma1, rho)
}





#' @title Calculates the first two moments of the bivariate distribution of NB_model and NB_all 
#' @param Y Vector of the binary response variable
#' @param pi Vector of predicted risks
#' @param z Decision threshold at which the NBs are calculated
#' @param weights Optinal - observation weights
#' @return Two means, two SDs, and one correlation coefficient. First element is for the model and second is for treating all
#' @export
calc_NB_moments <- function(Y,pi,z,weights=NULL){
  # set up
  n <- length(Y)
  a <- as.numeric(pi>z)
  if(is.null(weights))
  {
    rho <- mean(Y)
    TPR <- mean(Y*a)/rho
    FPR <- mean(a*(1-Y))/(1-rho)
  }
  else
  {
    rho <- stats::weighted.mean(Y,w = weights)
    TPR <- stats::weighted.mean(Y*a, w = weights)/rho
    FPR <- stats::weighted.mean(a*(1-Y), w=weights)/(1-rho) 
  }
  tz <- z/(1-z)
  
  # mean
  mu <- c(rho*TPR-tz*(1-rho)*FPR,rho-tz*(1-rho))
  
  # var
  sig_Y <- rho*(1-rho)/n
  sig_TPR <- rho*TPR*(1-rho*TPR)/n
  sig_FPR <- (1-rho)*FPR*(1-(1-rho)*FPR)/n
  
  sig11 <-  sig_TPR + tz^2 * sig_FPR + 2 * tz * rho * (1-rho) * TPR * FPR / n
  sig22 <- (1+tz)^2 * sig_Y
  sig12 <- sig_Y * (1+tz)* (TPR + tz * FPR)
  cor_coefficient <- sig12/(sqrt(sig11)*sqrt(sig22))
  
  # mean of NB_model, mean of NB_all,
  # var of NB_model, var of NB_all,
  # correlation of NB_model and NB_all
  return(c(mu1=mu[1],mu2=mu[2],sd1=sqrt(sig11),sd2=sqrt(sig22),rho=cor_coefficient))
}







