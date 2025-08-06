#' This function helps to check whether the models/patients are converging to
#' stable values. It should be used with the combined results
#' of the simulations/patients.

#' Inputs:

#' 1. results_combined from fun_results, either of patients or of simulations.

#' 2. Interval within which to check for convergence (5 models, 10 patients, etc.)
#' When the number of simulations/patients is divided by this number, it must lead to a whole number, not a fraction.
#' GOOD: 500 patients and interval of 100. 500 / 100 = 5 (whole number)
#' BAD: 500 patients and interval of 200. 500 / 200 = 2.5 (not whole number)

#' 3. Percentage within which cumulative mean should stay to count as converged.
#' The smaller this number, the stricter the criteria for convergence.

fun_diagnostics <- function(results_combined, interval, percentage){ 
  
  #===============================================================================
  # List
  #===============================================================================
  
  diagnostics <- vector(mode = "list", length = 2)
  names(diagnostics) <- c("cumulative_mean_interval","convergence")
  
  #===============================================================================
  # Cumulative Means
  #===============================================================================
  
  # For our combined results, we create an array to store the cumulative means
  # Row 1: Outcomes of simulation/patient 1
  # Row 2: Mean outcomes of simulations/patients 1 and 2
  # Row 3: Mean outcomes of simulations/patients 1 - 3
  # ...
  # Row 1'000: Mean outcomes of simulations/patients 1 - 1'000
  # Etc.
  
  cumulative_mean <- array(0,
                           dim=c(nrow(results_combined),ncol(results_combined)),
                           dimnames = list(c(1:nrow(results_combined)),
                                           names(results_combined)))
  
  for (i in 1:ncol(results_combined)) {
    cumulative_mean[,i] <- cumsum(results_combined[,i]) / seq_along(results_combined[,i])
  }
  
  # Convert cumulative means into data frame
  
  cumulative_mean <- data.frame(cumulative_mean)
  
  #===============================================================================
  # Convergence
  #===============================================================================
  
  # We define the intervals and the percentage change within which the outcomes count as converged
  
  X <- interval
  Y <- percentage
  
  intervals <- nrow(results_combined) / X  # This has to be a whole number
  
  # Create an array and store mean outcomes in intervals of X
  
  cumulative_mean_outcome_interval <- array(0,
                                            dim=c(intervals,ncol(results_combined)),
                                            dimnames = list(c(1:intervals),
                                                            names(results_combined)))
  
  
  for (j in 1:intervals) {
    for (i in 1:ncol(results_combined)) {
      cumulative_mean_outcome_interval[j,i] <- cumulative_mean[j*X,i]
    }
  }
  
  diagnostics$cumulative_mean_interval <- cumulative_mean_outcome_interval
  
  # Check whether mean outcome changed by more than Y% from last interval
  
  convergence_check <- array(0,
                             dim=c(intervals,ncol(results_combined)),
                             dimnames = list(c(1:intervals),
                                             names(results_combined)))
  
  # 1 => The mean outcome has changed by more than Y% if we also include the next interval in our calculation of the mean outcome
  # 0 => The mean outcome has converged (by our definition of a stable value)
  
  for (j in 2:intervals) {
    for (i in 1:ncol(results_combined)) {
      convergence_check[j,i] <- (abs(1-cumulative_mean_outcome_interval[j,i]/cumulative_mean_outcome_interval[j-1,i]) > Y)
    }
  }
  
  diagnostics$convergence <- convergence_check
  
  # We return our list of plots from the functions
  return(diagnostics)
  
}