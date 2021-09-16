

#' Title Estimate mean and variance of prediction based on model calibration output.
#'
#' @param calibVector Vector of predicted probability of risk per decile or percentile (e.g., from a calibration plot).
#'
#' @return Returns mean and variance of predictions based on the predicted probabilities.
#' @export
pred_summary_stat <- function(calibVector) {
  p <- mean(calibVector)
  v <- mean((calibVector - p) ^ 2)
  
  return(list(pred_mean = p, pred_var = v))
}




#' Title Update a prediction model for a binary outcome by multiplying a fixed odd-ratio to the predicted odds.
#'
#' @param p0 Mean of observed risk or predicted risk in development sample.
#' @param v Variance of predicted risk in development sample.
#' @param p1 Mean of observed risk in target population.
#'
#' @return Returns a correction factor that can be applied to the predicted 
#' odds in order to update the predictions for a new target population.
#' @export
odds_adjust <- function(p0, p1, v)
{

  if(v > p0 * (1 - p0)) stop("Variance cannot be larger than p0*(1-p0).")
  
  A <- p0 ^ 3 - p1 * p0 ^ 3
  B <- 3 * p1 * p0 ^ 3 - 2 * p0 ^ 3 - 3 * p1 * p0 ^ 2 + 2 * p0 ^ 2 - v
  C <- p0 ^ 3 - 3 * p1 * p0 ^ 3 + 6 * p1 * p0 ^ 2 - 2 * p0 ^ 2 + p0 - 3 * p1 * p0 + v
  D <- p1 * p0 ^ 3 - 3 * p1 * p0 ^ 2 + 3 * p1 * p0 - p1
  
  res <- cubic(c(A, B, C, D))
  res <- Re(res[which(Im(res) == 0)]) #Remove non-reals
  res <- res[which(res > 0)] #Remove negatives
  res <- res[which(sign(log(res)) == sign(log(p1 / p0)))] #Removes ORs in the wrong direction
  res
}


