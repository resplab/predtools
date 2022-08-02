#' @title Calculates the expected value of the maximum of two random variables with zero-truncated bivariate normal distirbution
#' Takes a vector of mean and a 2X2 covariance matrix
#' @param mu1 Mean of the first distribution
#' @param mu2 Mean of the second distribution
#' @param sd1 SD of the first distribution
#' @param sd2 SD of the second distribution
#' @param rho Correlation coefficient of the two random variables
#' @return A scalar value for the expected value
#' @export

mu_max_trunc_bvn <-
  function(mu1, mu2, sd1, sd2, rho) {
    
    f1 <-  function(mu1, mu2, sd1, sd2, rho) {
      tmp1 <- sd1 - rho * sd2
      tmp2 <-
        (-sd1 * mu2 + rho * sd2 * mu1) / (sd1 * sd2 * sqrt(1 - rho ^
                                                                 2))
      mu1 * (as.numeric(tmp1 > 0) +  0 * as.numeric(tmp1 == 0) * pnorm(tmp2)) -
        pnorm(tmp2) * (-sd1 * dnorm(-mu1 / sd1) + mu1 * pnorm(-mu1 / sd1))
    }
    
    f2 <- function(mu1, mu2, sd1, sd2, rho) {
      tmp1 <- (sd1 - rho * sd2)
      alpha_num <- (sd1 * mu2 - rho * sd2 * mu1)
      beta_num <- sd1 * sd2 * sqrt(1 - rho ^ 2)
      alpha <- alpha_num / tmp1
      beta <- beta_num / tmp1
      
      if (tmp1 > 0) {
        a <- (tmp1 * mu1 - alpha_num) / beta_num
        b <- tmp1 * sd1 / beta_num
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
        a <- (alpha - mu1) / sd1
        b <- beta / sd1
        T2_t <- sqrt(1 + b ^ 2)
        
        if (tmp1 < 1e-10) {
          T2 <-
            -sd1 / b * dnorm(alpha_num / beta_num) * (1 - pnorm(-mu1 / sd1))
        } else{
          T2 <-
            -sd1 / T2_t * dnorm(a / T2_t) *
            (1 - pnorm(-T2_t * alpha_num / beta_num + a * b / T2_t))
        }
        
        return(T1 + T2)
      }
      else{
        beta <- abs(beta)
        
        a <- (abs(tmp1) * mu1 + alpha_num) / beta_num
        b <- abs(tmp1) * sd1 / beta_num
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
        
        a <- (alpha - mu1) / sd1
        b <- beta / sd1
        T2_t <- sqrt(1 + b ^ 2)
        
        if (abs(tmp1) < 1e-10) {
          T2 <-
            sd1 / b * dnorm(-alpha_num / beta_num) * (1 - pnorm(-mu1 / sd1))
        } else{
          T2 <-
            sd1 / T2_t * dnorm(a / T2_t) * (1 - pnorm(T2_t * alpha_num / beta_num +
                                                         a * b / T2_t))
        }
        
        return(T1 + T2)
        
      }
    }
    
    f1(mu1, mu2, sd1, sd2, rho) + f1(mu2, mu1, sd2, sd1, rho) -
      f2(mu1, mu2, sd1, sd2, rho) - f2(mu2, mu1, sd2, sd1, rho)
  }
