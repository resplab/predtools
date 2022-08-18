evpi_val <- function(Y, pi, method=c("bootstrap","bayesian bootstrap","asymptotic"), n_sim=1000, zs=(0:99)/100, weights=NULL)
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
          {ENB_perfect[j] <- do.call(mu_max_trunc_bvn,as.list(parms))}
          , error=function(cond) {
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
    rho <- weighted.mean(Y,w = weights)
    TPR <- weighted.mean(Y*a, w = weights)/rho
    FPR <- weighted.mean(a*(1-Y), w=weights)/(1-rho) 
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







