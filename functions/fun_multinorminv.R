#' Function to randomly draw coefficients based on Cholesky decomposition of covariance matrix from Ara and Brazier (2010)
#' Function Input format:
#' Mu: Mu must be a vector of 4 values, the coefficients from the Ara and Brazier regression
#' S: S must be a 4x4 matrix, with the covariance from the Ara and Brazier regression
#' Z: Z must be a vector of 4 random values between 0 and 1

fun_multinorminv <- function(Mu, S, Z) {
  
  # We use the same notation as MULTINORMINV function
  # Mu: Coefficients from regression model
  # S: Covariance Matrix
  # L: Cholesky decomposition of covariance matrix
  # Z: Standardized random numbers based on normal distribution
  # W: Randomly drawn coefficients from multivariate normal distribution
  
  L <- chol(S)
  
  W <- Mu + L%*% Z
  
  return(W)
}
