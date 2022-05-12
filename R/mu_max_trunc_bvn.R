#' @title Calculates the expected value of the maximum of two random variables with zero-truncated bivariate normal distirbution
#' Takes a vector of mean and a 2X2 covariance matrix
#' @param mu A vector of means (of length 2)
#' @param sigma A 2X2 covariance matrix.
#' @return A scalar value for the expected value
#' @export

mu_max_trunc_bvn <-
  function(mu, sigma) {
    mu1 <- mu[1]
    mu2 <- mu[2]
    sig1 <- sqrt(sigma[1,1])
    sig2 <- sqrt(sigma[2,2])
    rho <- sigma[1,2]/sqrt(sig1*sig2)
    
    f1 <-  function(mu1, mu2, sig1, sig2, rho) {
      tmp1 <- sig1 - rho * sig2
      tmp2 <-
        (-sig1 * mu2 + rho * sig2 * mu1) / (sig1 * sig2 * sqrt(1 - rho ^
                                                                 2))
      mu1 * (as.numeric(tmp1 > 0) +  0 * as.numeric(tmp1 == 0) * pnorm(tmp2)) -
        pnorm(tmp2) * (-sig1 * dnorm(-mu1 / sig1) + mu1 * pnorm(-mu1 / sig1))
    }
    
    f2 <- function(mu1, mu2, sig1, sig2, rho) {
      tmp1 <- (sig1 - rho * sig2)
      alpha_num <- (sig1 * mu2 - rho * sig2 * mu1)
      beta_num <- sig1 * sig2 * sqrt(1 - rho ^ 2)
      alpha <- alpha_num / tmp1
      beta <- beta_num / tmp1
      
      if (tmp1 > 0) {
        a <- (tmp1 * mu1 - alpha_num) / beta_num
        b <- tmp1 * sig1 / beta_num
        T1_rho <- -1 / sqrt(1 + b ^ 2)
        
        if (tmp1 < 1e-10) {
          T1_1 <- pnorm(alpha_num / beta_num)
        } else{
          T1_1 <-  pnorm((-a / b) / sqrt(1 + (1 / b) ^ 2))
        }
        
        T1 <-
          mu1 * (T1_1  - as.numeric(pmvnorm(
            lower = c(-Inf,-Inf),
            upper = c(-a / sqrt(1 + b ^ 2),-alpha_num / beta_num),
            mean =
              c(0, 0),
            sigma = matrix(c(1, T1_rho, T1_rho, 1), 2)
          )[[1]]))
        a <- (alpha - mu1) / sig1
        b <- beta / sig1
        T2_t <- sqrt(1 + b ^ 2)
        
        if (tmp1 < 1e-10) {
          T2 <-
            -sig1 / b * dnorm(alpha_num / beta_num) * (1 - pnorm(-mu1 / sig1))
        } else{
          T2 <-
            -sig1 / T2_t * dnorm(a / T2_t) *
            (1 - pnorm(-T2_t * alpha_num / beta_num + a * b / T2_t))
        }
        
        return(T1 + T2)
      }
      else{
        beta <- abs(beta)
        
        a <- (abs(tmp1) * mu1 + alpha_num) / beta_num
        b <- abs(tmp1) * sig1 / beta_num
        T1_rho <- -1 / sqrt(1 + b ^ 2)
        
        if (abs(tmp1) < 1e-10) {
          T1_1 <- pnorm(-alpha_num / beta_num)
        } else{
          T1_1 <-  pnorm((-a / b) / sqrt(1 + (1 / b) ^ 2))
        }
        
        T1 <-
          -mu1 * (T1_1  - as.numeric(pmvnorm(
            lower = c(-Inf,-Inf),
            upper = c(-a / sqrt(1 + b ^ 2), alpha_num / beta_num),
            mean =
              c(0, 0),
            sigma = matrix(c(1, T1_rho, T1_rho, 1), 2)
          )[[1]]))
        
        a <- (alpha - mu1) / sig1
        b <- beta / sig1
        T2_t <- sqrt(1 + b ^ 2)
        
        if (abs(tmp1) < 1e-10) {
          T2 <-
            sig1 / b * dnorm(-alpha_num / beta_num) * (1 - pnorm(-mu1 / sig1))
        } else{
          T2 <-
            sig1 / T2_t * dnorm(a / T2_t) * (1 - pnorm(T2_t * alpha_num / beta_num +
                                                         a * b / T2_t))
        }
        
        return(T1 + T2)
        
      }
    }
    
    f1(mu1, mu2, sig1, sig2, rho) + f1(mu2, mu1, sig2, sig1, rho) -
      f2(mu1, mu2, sig1, sig2, rho) - f2(mu2, mu1, sig2, sig1, rho)
  }