################################################################################
#                                                                              #
# Study: Cost-Effectiveness Analysis of Gene Therapy for Haemophilia B         #
# Design: Cost-Effectiveness Model using Bleeds as a continuous variable       #
# Outcome: Costs, QALYS, ICER                                                  #
# Task: Univariate sensitivity analysis                                        #
# Author: Niklaus Meier                                                        #
# R version: 4.2.1                                                             #
#                                                                              #
################################################################################

start.time <- Sys.time()

################################################################################
#                                                                              #
# Setup                                                                        #
#                                                                              #
################################################################################

## Load data
load(file = paste0(directories$dir_dat_deriv, "/Germany_life_tables_weight.RData"))

################################################################################
#                                                                              #
# Baseline Analysis                                                            #
#                                                                              #
################################################################################

# We run a standard deterministic model to set our baseline

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = 10000,
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)
model <- fun_gen_pop(model = model, mode = "heterogeneous")

# We save our generated population so we can re-use it later
population <- model$sim

# We also save our baseline parameters
param_sample <- model$params$prob_sample

# We run the model as normal
model <- fun_ETRANACOGENE_abr(model = model, mode = "random")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)

# We save the results of our baseline model
baseline_patient_results <- fun_results(model, mode = "patient")
baseline_model_results   <- fun_results(model, mode = "simulation")

################################################################################
#                                                                              #
# Sensitivity Parameters and Ranges                                            #
#                                                                              #
################################################################################

sens_analysis <- list()

#===============================================================================
# Parameters
#===============================================================================

# We select and name the parameters for which we want to perform our sensitivity analysis
# These names NEED to correspond exactly to those in fun_param_sample

# We only include parameters that are probabilistic and not related to population generation

#-------------------------------------------------------------------------------
# Resource Use
#-------------------------------------------------------------------------------

sens_analysis$res_hosp_prob$name <- "res_hosp_prob"
sens_analysis$res_hosp_LOS$name  <- "res_hosp_LOS"

#-------------------------------------------------------------------------------
# Utilities
#-------------------------------------------------------------------------------

sens_analysis$util_baseline$name     <- "util_baseline"
sens_analysis$util_male$name         <- "util_male"
sens_analysis$util_age$name          <- "util_age"
sens_analysis$util_agesq$name        <- "util_agesq"
sens_analysis$util_PS_13_to_21$name  <- "util_PS_13_to_21"
sens_analysis$util_PS_22_plus$name   <- "util_PS_22_plus"
sens_analysis$util_joint_bleed$name  <- "util_joint_bleed"
sens_analysis$util_other_bleed$name  <- "util_other_bleed"
sens_analysis$util_infusion$name     <- "util_infusion"
sens_analysis$util_surgery$name      <- "util_surgery"

#-------------------------------------------------------------------------------
# Clinical Parameters
#-------------------------------------------------------------------------------

sens_analysis$share_joint_bleeds$name    <- "share_joint_bleeds"
sens_analysis$bleed_mortality_ratio$name <- "bleed_mortality_ratio"
sens_analysis$PS_bleeds$name             <- "PS_bleeds"

#-------------------------------------------------------------------------------
# ETRANACOGENE
#-------------------------------------------------------------------------------

sens_analysis$ETRANACOGENE_relative_bleed_reduction$name            <- "ETRANACOGENE_relative_bleed_reduction"
sens_analysis$ETRANACOGENE_max_bleed_reduction_duration$name        <- "ETRANACOGENE_max_bleed_reduction_duration"
sens_analysis$ETRANACOGENE_bleed_increase_per_year$name             <- "ETRANACOGENE_bleed_increase_per_year"
sens_analysis$ETRANACOGENE_success_prob$name                        <- "ETRANACOGENE_success_prob"

#-------------------------------------------------------------------------------
# Prophylaxis
#-------------------------------------------------------------------------------

sens_analysis$PROPHYLAXIS_relative_bleed_reduction$name          <- "PROPHYLAXIS_relative_bleed_reduction"

#===============================================================================
# Value Ranges
#===============================================================================

# We define low and high bounds for our parameters, to see how outcomes change based on this.

# For those values based on distribution, we choose the upper and lower bound based on the
# 2.5% and 97.5% quantile of their distribution, based on the parameters of the probabilistic distribution.

# For those values based on assumptions and uniform distributions, we choose the upper and lower bound based on the
# lowest and highest values in our PSA.

#-------------------------------------------------------------------------------
# Resource Use
#-------------------------------------------------------------------------------

sens_analysis$res_hosp_prob$range <- c(qbeta(0.025, shape1 =  model$params$resources$hosp_ratio$hosps + 1, shape2 = model$params$resources$hosp_ratio$bleeds - model$params$resources$hosp_ratio$hosps + 1),
                                       qbeta(0.975, shape1 =  model$params$resources$hosp_ratio$hosps + 1, shape2 = model$params$resources$hosp_ratio$bleeds - model$params$resources$hosp_ratio$hosps + 1))

sens_analysis$res_hosp_LOS$range  <- c(qgamma(0.025, shape = model$params$resources$hosp_LOS$alpha , scale = model$params$resources$hosp_LOS$beta),
                                       qgamma(0.975, shape = model$params$resources$hosp_LOS$alpha, scale =  model$params$resources$hosp_LOS$beta))

#-------------------------------------------------------------------------------
# Utilities
#-------------------------------------------------------------------------------

# Our utilities for baseline, male, age, and agesq are all mutually dependent on each other due to the covariance matrix.
# Therefore, we can't sample them independently.
# to get the 2.5% and 97.% quantile, we sample them 10'000 times and determine the quantiles based on this

temp_utils <- matrix(nrow = 4, ncol = 10000)

for (i in 1:ncol(temp_utils)) {
  temp_utils[,i] <- fun_multinorminv(model$params$utilities$mu,
                                     model$params$utilities$cov_matrix,
                                     qnorm(runif(4), mean = 0, sd = 1))
}

# We assign the these defined quantiles directly as the lower and upper bound for these utilities

sens_analysis$util_male$range         <- as.vector(quantile(temp_utils[1,], probs =  c(0.025, 0.975)))

sens_analysis$util_age$range          <- as.vector(quantile(temp_utils[2,], probs =  c(0.025, 0.975)))

sens_analysis$util_agesq$range        <- as.vector(quantile(temp_utils[3,], probs =  c(0.025, 0.975)))

sens_analysis$util_baseline$range     <- as.vector(quantile(temp_utils[4,], probs =  c(0.025, 0.975)))

rm(temp_utils)

# For the other utilities we can find the quantiles directly from the distributions

sens_analysis$util_PS_13_to_21$range  <- c(-qgamma(0.025, shape = model$params$utilities$PS_13_to_21$alpha, scale = model$params$utilities$PS_13_to_21$beta),
                                           -qgamma(0.975, shape = model$params$utilities$PS_13_to_21$alpha, scale = model$params$utilities$PS_13_to_21$beta))

sens_analysis$util_PS_22_plus$range   <- c(-qgamma(0.025, shape = model$params$utilities$PS_22_plus$alpha, scale = model$params$utilities$PS_22_plus$beta),
                                           -qgamma(0.975, shape = model$params$utilities$PS_22_plus$alpha, scale = model$params$utilities$PS_22_plus$beta))

sens_analysis$util_joint_bleed$range  <- c(-qgamma(0.025, shape = model$params$utilities$bleed$alpha, scale = model$params$utilities$bleed$beta),
                                           -qgamma(0.975, shape = model$params$utilities$bleed$alpha, scale = model$params$utilities$bleed$beta))

sens_analysis$util_other_bleed$range  <- c(-qgamma(0.025, shape = model$params$utilities$bleed$alpha, scale = model$params$utilities$bleed$beta),
                                           -qgamma(0.975, shape = model$params$utilities$bleed$alpha, scale = model$params$utilities$bleed$beta))

sens_analysis$util_infusion$range     <- c(-qgamma(0.025, shape = model$params$utilities$infusion$alpha, scale = model$params$utilities$infusion$beta),
                                           -qgamma(0.975, shape = model$params$utilities$infusion$alpha, scale = model$params$utilities$infusion$beta))

sens_analysis$util_surgery$range      <- c(-qgamma(0.025, shape = model$params$utilities$surgery$alpha, scale = model$params$utilities$surgery$beta) * 27/model$struct$cyc2day,
                                           -qgamma(0.975, shape = model$params$utilities$surgery$alpha, scale = model$params$utilities$surgery$beta) * 27/model$struct$cyc2day)

#-------------------------------------------------------------------------------
# Clinical Parameters
#-------------------------------------------------------------------------------

sens_analysis$share_joint_bleeds$range    <- c(qbeta(0.025, shape1 = model$params$joint_bleeds$events + 1, shape2 =  model$params$joint_bleeds$nsample - model$params$joint_bleeds$events + 1),
                                               qbeta(0.975, shape1 = model$params$joint_bleeds$events + 1, shape2 =  model$params$joint_bleeds$nsample - model$params$joint_bleeds$events + 1))

# Share other bleeds is directly determined by joint bleeds and not sampled itself. Since we can't sample both at the same time for our univariate sensitivity analysis,
# we make this correction down in the loop below.

sens_analysis$bleed_mortality_ratio$range <- c(qgamma(0.025, shape = model$params$bleed_mort$alpha, scale = model$params$bleed_mort$beta),
                                               qgamma(0.975, shape = model$params$bleed_mort$alpha, scale = model$params$bleed_mort$beta))

sens_analysis$PS_bleeds$range             <- c(qgamma(0.025, shape = model$params$pettersson$alpha, scale = model$params$pettersson$beta),
                                               qgamma(0.975, shape = model$params$pettersson$alpha, scale = model$params$pettersson$beta))

#-------------------------------------------------------------------------------
# ETRANACOGENE
#-------------------------------------------------------------------------------

# Take quantiles of ABR and calculate bleed reduction based on ABR quantiles

sens_analysis$ETRANACOGENE_relative_bleed_reduction$range           <- 1 - (c(qgamma(0.025, shape = model$params$ETRANACOGENE$abr$alpha, scale = model$params$ETRANACOGENE$abr$beta),
                                                                           qgamma(0.975, shape = model$params$ETRANACOGENE$abr$alpha, scale = model$params$ETRANACOGENE$abr$beta))
                                                                         /model$params$baseline_abr$mean)

sens_analysis$ETRANACOGENE_max_bleed_reduction_duration$range       <- c(model$params$ETRANACOGENE$max_bleed_reduction_duration$l_limit, 
                                                                      model$params$ETRANACOGENE$max_bleed_reduction_duration$u_limit)

sens_analysis$ETRANACOGENE_bleed_increase_per_year$range            <- c(model$params$ETRANACOGENE$bleed_increase_per_year$l_limit, 
                                                                      model$params$ETRANACOGENE$bleed_increase_per_year$u_limit)

sens_analysis$ETRANACOGENE_success_prob$range                       <- c(qbeta(0.025, shape1 = model$params$ETRANACOGENE$success_prob$events + 1, shape2 = model$params$ETRANACOGENE$success_prob$nsample - model$params$ETRANACOGENE$success_prob$events + 1), 
                                                                      qbeta(0.975, shape1 = model$params$ETRANACOGENE$success_prob$events + 1, shape2 = model$params$ETRANACOGENE$success_prob$nsample - model$params$ETRANACOGENE$success_prob$events + 1))

#-------------------------------------------------------------------------------
# Prophylaxis
#-------------------------------------------------------------------------------

# Take quantiles of ABR and calculate bleed reduction based on ABR quantiles

sens_analysis$PROPHYLAXIS_relative_bleed_reduction$range          <- 1 - (c(qgamma(0.025, shape = model$params$PROPHYLAXIS$abr$alpha, scale = model$params$PROPHYLAXIS$abr$beta),
                                                                            qgamma(0.975, shape = model$params$PROPHYLAXIS$abr$alpha, scale = model$params$PROPHYLAXIS$abr$beta))
                                                                          /model$params$baseline_abr$mean)

################################################################################
#                                                                              #
# Sensitivity Analysis                                                         #
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Loop
#-------------------------------------------------------------------------------

#' The univariate sensitivity analysis can require a lot of memory if many 
#' simulations are generated with a large number of patients.
#' Therefore, it may be necessary to run some simulations, save the outputs,
#' and then delete the model from the environment before running more simulations.
#' We do this by running multiple loops with 1 simulation per loop.

## We prepare an empty data frame to save our results across multiple loops
# of our probabilistic sensitivity analysis

gc()
model_results_combined   <- data.frame()
patient_results_combined <- data.frame()

for (j in 1:length(sens_analysis)) {
  
  print(paste0("Parameter ", j, "/", length(sens_analysis), ": " ,  sens_analysis[[j]][["name"]]))
  
  for (i in 1:length(sens_analysis[[j]][["range"]])) {
    
    # We reset the population in our model to that of our baseline model
    
    model$sim <- population
    
    # We reset parameters in model
    
    model$params$prob_sample <- param_sample
    
    # We assign the correct parameter value for our sensitivity analysis
    
    model[["params"]][["prob_sample"]][,paste0(sens_analysis[[j]][["name"]])] <- sens_analysis[[j]][["range"]][[i]]
    
    # Adjust other bleeds which are directly dependent on joint bleeds
    
    model[["params"]][["prob_sample"]][,"share_other_bleeds"] <- 1 - model[["params"]][["prob_sample"]][,"share_joint_bleeds"]
    
    # Run Model
    
    model <- fun_ETRANACOGENE_abr(model = model, mode = "random")
    model <- fun_PROPHYLAXIS_abr(model = model)
    model <- fun_bleeds_death(model = model, mode = "expected_value")
    model <- fun_resources(model = model)
    model <- fun_utility(model = model)
    model <- fun_costs(model = model)
    
    ################################################################################
    #                                                                              #
    # Saving Results                                                               #
    #                                                                              #
    ################################################################################
    
    # Save results of model to data frame
    patient_results <- fun_results(model, mode = "patient")
    model_results   <- fun_results(model, mode = "simulation")
    
    patient_results[,"sens_analysis_param"] <- sens_analysis[[j]][["name"]]
    model_results[,"sens_analysis_param"]  <- sens_analysis[[j]][["name"]]
    
    # Combine results of this loop with results of previous loops
    patient_results_combined <- rbind(patient_results_combined, patient_results)
    model_results_combined   <- rbind(model_results_combined, model_results)
    
    # Print progress
    print(paste0(i, "/", length(sens_analysis[[j]][["range"]]), " simulations complete"))
    
    # Clean up memory
    gc()
  }
  
}

# We assign a number to each patient.

patient_results_combined[,"patient_number"] <- 1:nrow(patient_results_combined)

# Save sensitivity analysis results to list based on parameter

for (j in 1:length(sens_analysis)) {
  sens_analysis[[j]][["patient_results"]] <- subset(patient_results_combined,sens_analysis_param == paste0(sens_analysis[[j]][["name"]]))
  sens_analysis[[j]][["model_results"]]   <- subset(model_results_combined,sens_analysis_param == paste0(sens_analysis[[j]][["name"]]))
}

################################################################################
#                                                                              #
# Save Results                                                                 #
#                                                                              #
################################################################################

save(patient_results_combined, file = paste0(directories$dir_dat_deriv, '/sens_analysis_patients.Rdata'))
save(model_results_combined, file = paste0(directories$dir_dat_deriv, '/sens_analysis_model.Rdata'))
save(sens_analysis, file = paste0(directories$dir_dat_deriv, '/sens_analysis_inputs_results.Rdata'))
save(baseline_model_results, file = paste0(directories$dir_dat_deriv, '/sens_analysis_baseline.Rdata'))
save(population, file = paste0(directories$dir_dat_deriv, '/sens_analysis_population.Rdata'))

################################################################################
#                                                                              #
# Cleanup                                                                      #
#                                                                              #
################################################################################

rm(model, i, life_tables_weight,
   baseline_model_results, baseline_patient_results, 
   model_results, patient_results,
   param_sample, population,
   patient_results_combined, model_results_combined, sens_analysis)

################################################################################
#                                                                              #
# REPORT RUNTIME                                                               #
#                                                                              #
################################################################################

end.time <- Sys.time()

print(paste('Run time =', round(as.numeric(end.time, units = "secs") - as.numeric(start.time, units = "secs"), 2)/60, 'minutes', sep = ' '))

gc()