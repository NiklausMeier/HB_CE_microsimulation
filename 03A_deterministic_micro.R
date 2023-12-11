################################################################################
#                                                                              #
# Study: Cost-Effectiveness Analysis of Gene Therapy for Haemophilia B         #
# Design: Cost-Effectiveness Model using Bleeds as a continuous variable       #
# Outcome: Costs, QALYS, ICER                                                  #
# Task: Run deterministic version of model with baseline parameters            #
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
# Prepare Model                                                                #
#                                                                              #
################################################################################

#' To prepare the model and set parameters, the next functions must be applied in the following order:
#' 1. Generate model object: fun_model_setup
#' 2. Set baseline parameters: fun_parameter_baseline
#' 3. Sample probabilistic parameters: fun_param_sample
#' 4. Set scenarios: fun_scenarios

#===============================================================================
# Generate model object
#===============================================================================

## Choose cycle length, number of cycles, discount rate and the number of patients
# Cycle length should be no shorter than 27 days, due to the calculation of disutilities for some events

# We also need to choose whether the model is "deterministic" or "probabilistic"
# This determines how parameters are sampled when using the function "fun_param_sample"

# As treatments we include ETRANACOGENE and PROPHYLAXIS
# The order determines which number the treatments receive in the data array, and should therefore not be changed

# The threshold determines how many units of a currency [EUR] we are willing to pay per QALY
# This is used to determine cost-effectiveness and calculate the Net Monetary Benefit (NMB)

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = 10000,
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

#===============================================================================
# Set baseline parameters
#===============================================================================

# We use this function to set baseline parameters. 
# This includes means, but also measures of uncertainty (variance, standard deviation, or standard error).
# These values are used for the deterministic as well as the probabilistic model.
# To change a baseline value, the function should be edited directly, to set a new baseline.
# To conduct e.g. sensitivity analysis, individual values can then be edited directly after running the function.

model <- fun_parameter_baseline(model = model)

#===============================================================================
# Sample probabilistic parameters
#===============================================================================

model <- fun_param_sample(model = model)

#===============================================================================
# Scenarios
#===============================================================================

model <- fun_scenarios(model = model)

################################################################################
#                                                                              #
# Run Model                                                                    #
#                                                                              #
################################################################################

#' Using the model structure, the next functions must be applied in the following order:
#' 1. Generate population: fun_gen_pop
#' 2. Determine ABR based on treatment: fun_ETRANACOGENE_abr, or fun_PROPHYLAXIS_abr
#' 3. Determine bleeding and death: fun_bleeds_death
#' 4. Determine resource use and costs: fun_resources
#' 5. Determine utilities and QALYs: fun_utility

#===============================================================================
# Generate Population
#===============================================================================

model <- fun_gen_pop(model = model, mode = "heterogeneous")

# We save our generated population so we can re-use it later
population <- model$sim

#===============================================================================
# Treatment
#===============================================================================

model <- fun_ETRANACOGENE_abr(model = model, mode = "random")
model <- fun_PROPHYLAXIS_abr(model = model)

#===============================================================================
# Bleeding and Death
#===============================================================================

model <- fun_bleeds_death(model = model, mode = "expected_value")

#===============================================================================
# Resources
#===============================================================================

model <- fun_resources(model = model)

#===============================================================================
# Utility
#===============================================================================

model <- fun_utility(model = model)

#===============================================================================
# Costs
#===============================================================================

model <- fun_costs(model = model)

#===============================================================================
# Results
#===============================================================================

patient_results <- fun_results(model, mode = "patient")
model_results   <- fun_results(model, mode = "simulation")

#===============================================================================
# Diagnostics
#===============================================================================

patient_diag <- fun_diagnostics(patient_results, interval = 1000, percentage = 0.01)

#===============================================================================
# Outcomes over time
#===============================================================================

# Outcomes we want to track over time
outcomes_over_time <- c("Years",
                        
                        "ETRANACOGENE_LYs_cycle_undisc",   "PROPHYLAXIS_LYs_cycle_undisc",
                        "ETRANACOGENE_QALYs_cycle_undisc", "PROPHYLAXIS_QALYs_cycle_undisc",
                        "ETRANACOGENE_cost_cycle_undisc",  "PROPHYLAXIS_cost_cycle_undisc",
                        
                        "ETRANACOGENE_LYs_sum_undisc",   "PROPHYLAXIS_LYs_sum_undisc",
                        "ETRANACOGENE_QALYs_sum_undisc", "PROPHYLAXIS_QALYs_sum_undisc",
                        "ETRANACOGENE_cost_sum_undisc",  "PROPHYLAXIS_cost_sum_undisc",
                        
                        "disc_factor",
                        
                        "ETRANACOGENE_LYs_cycle_disc",   "PROPHYLAXIS_LYs_cycle_disc",
                        "ETRANACOGENE_QALYs_cycle_disc", "PROPHYLAXIS_QALYs_cycle_disc",
                        "ETRANACOGENE_cost_cycle_disc",  "PROPHYLAXIS_cost_cycle_disc",
                        
                        "ETRANACOGENE_LYs_sum_disc",   "PROPHYLAXIS_LYs_sum_disc",
                        "ETRANACOGENE_QALYs_sum_disc", "PROPHYLAXIS_QALYs_sum_disc",
                        "ETRANACOGENE_cost_sum_disc",  "PROPHYLAXIS_cost_sum_disc"
                        
)

# Create array with one row per cycle
outcomes_over_time_array <- setNames(data.frame(matrix(data = 0,
                                                       nrow = model$struct$ncycle,
                                                       ncol = length(outcomes_over_time))), 
                                     outcomes_over_time)

# Discount factors
outcomes_over_time_array[,"disc_factor"] <- model[["disc_factors"]][["disc_factors"]]


for (i in 1:model$struct$ncycle){
  
  # Number of years passed at end of cycle
  
  outcomes_over_time_array[i,"Years"]                               <- i * model$struct$cyc2day / settings$time$yr2day
  
  # Undiscounted outcomes, per cycle, and as a cumulative sum
  
  outcomes_over_time_array[i,"ETRANACOGENE_LYs_cycle_undisc"]       <- mean(model[["sim"]][["ETRANACOGENE"]][[1]][i,"LYs_undisc",])
  outcomes_over_time_array[i,"ETRANACOGENE_LYs_sum_undisc"]         <- sum(outcomes_over_time_array[1:i,"ETRANACOGENE_LYs_cycle_undisc"])
  
  outcomes_over_time_array[i,"PROPHYLAXIS_LYs_cycle_undisc"]        <- mean(model[["sim"]][["PROPHYLAXIS"]][[1]][i,"LYs_undisc",])
  outcomes_over_time_array[i,"PROPHYLAXIS_LYs_sum_undisc"]          <- sum(outcomes_over_time_array[1:i,"PROPHYLAXIS_LYs_cycle_undisc"])
  
  outcomes_over_time_array[i,"ETRANACOGENE_QALYs_cycle_undisc"]     <- mean(model[["sim"]][["ETRANACOGENE"]][[1]][i,"QALYs_undisc",])
  outcomes_over_time_array[i,"ETRANACOGENE_QALYs_sum_undisc"]       <- sum(outcomes_over_time_array[1:i,"ETRANACOGENE_QALYs_cycle_undisc"])
  
  outcomes_over_time_array[i,"PROPHYLAXIS_QALYs_cycle_undisc"]      <- mean(model[["sim"]][["PROPHYLAXIS"]][[1]][i,"QALYs_undisc",])
  outcomes_over_time_array[i,"PROPHYLAXIS_QALYs_sum_undisc"]        <- sum(outcomes_over_time_array[1:i,"PROPHYLAXIS_QALYs_cycle_undisc"])
  
  outcomes_over_time_array[i,"ETRANACOGENE_cost_cycle_undisc"]      <- mean(model[["sim"]][["ETRANACOGENE"]][[1]][i,"cost_undisc",])
  outcomes_over_time_array[i,"ETRANACOGENE_cost_sum_undisc"]        <- sum(outcomes_over_time_array[1:i,"ETRANACOGENE_cost_cycle_undisc"])
  
  outcomes_over_time_array[i,"PROPHYLAXIS_cost_cycle_undisc"]       <- mean(model[["sim"]][["PROPHYLAXIS"]][[1]][i,"cost_undisc",])
  outcomes_over_time_array[i,"PROPHYLAXIS_cost_sum_undisc"]         <- sum(outcomes_over_time_array[1:i,"PROPHYLAXIS_cost_cycle_undisc"])
  
  # Discounted outcomes, per cycle, and as a cumulative sum
  
  outcomes_over_time_array[i,"ETRANACOGENE_LYs_cycle_disc"]       <- outcomes_over_time_array[i,"ETRANACOGENE_LYs_cycle_undisc"]   * outcomes_over_time_array[i,"disc_factor"]
  outcomes_over_time_array[i,"ETRANACOGENE_LYs_sum_disc"]         <- sum(outcomes_over_time_array[1:i,"ETRANACOGENE_LYs_cycle_disc"])
  
  outcomes_over_time_array[i,"PROPHYLAXIS_LYs_cycle_disc"]        <- outcomes_over_time_array[i,"PROPHYLAXIS_LYs_cycle_undisc"]    * outcomes_over_time_array[i,"disc_factor"]
  outcomes_over_time_array[i,"PROPHYLAXIS_LYs_sum_disc"]          <- sum(outcomes_over_time_array[1:i,"PROPHYLAXIS_LYs_cycle_disc"])
  
  outcomes_over_time_array[i,"ETRANACOGENE_QALYs_cycle_disc"]     <- outcomes_over_time_array[i,"ETRANACOGENE_QALYs_cycle_undisc"] * outcomes_over_time_array[i,"disc_factor"]
  outcomes_over_time_array[i,"ETRANACOGENE_QALYs_sum_disc"]       <- sum(outcomes_over_time_array[1:i,"ETRANACOGENE_QALYs_cycle_disc"])
  
  outcomes_over_time_array[i,"PROPHYLAXIS_QALYs_cycle_disc"]      <- outcomes_over_time_array[i,"PROPHYLAXIS_QALYs_cycle_undisc"]  * outcomes_over_time_array[i,"disc_factor"]
  outcomes_over_time_array[i,"PROPHYLAXIS_QALYs_sum_disc"]        <- sum(outcomes_over_time_array[1:i,"PROPHYLAXIS_QALYs_cycle_disc"])
  
  outcomes_over_time_array[i,"ETRANACOGENE_cost_cycle_disc"]      <- outcomes_over_time_array[i,"ETRANACOGENE_cost_cycle_undisc"]  * outcomes_over_time_array[i,"disc_factor"]
  outcomes_over_time_array[i,"ETRANACOGENE_cost_sum_disc"]        <- sum(outcomes_over_time_array[1:i,"ETRANACOGENE_cost_cycle_disc"])
  
  outcomes_over_time_array[i,"PROPHYLAXIS_cost_cycle_disc"]       <- outcomes_over_time_array[i,"PROPHYLAXIS_cost_cycle_undisc"]   * outcomes_over_time_array[i,"disc_factor"]
  outcomes_over_time_array[i,"PROPHYLAXIS_cost_sum_disc"]         <- sum(outcomes_over_time_array[1:i,"PROPHYLAXIS_cost_cycle_disc"])
  
}

################################################################################
#                                                                              #
# Code to double-check generated population and model                          #
#                                                                              #
################################################################################

## Set i and j to 1
# i <- 1
# j <- 1

## Check all generated values of a single patient
# model[["sim"]][[j]][[i]][,,1]

## Check ... of all patients in the first cycle in a simulation
# model[["sim"]][[j]][[i]][1,"...",]

## Check ... of a single patient for all cycle in a simulation
# model[["sim"]][[j]][[i]][,"...",1]

##  Check starting ... of all patients in a simulation, and see if they are the same in treatment arms
# model[["sim"]][[j]][[i]][1,"...",]
# model[["sim"]][[j]][[i]][1,"...",] == model[["sim"]][[j+1]][[i]][1,"...",]

## Create data table for a patient in simulation ...
# table <- data.table(model[["sim"]][[j]][[1]][,,...])

################################################################################
#                                                                              #
# Save Results                                                                 #
#                                                                              #
################################################################################

save(model, file = paste0(directories$dir_dat_deriv, '/deterministic_model.Rdata'))
save(patient_results, file = paste0(directories$dir_dat_deriv, '/deterministic_patient_results.Rdata'))
save(model_results, file = paste0(directories$dir_dat_deriv, '/deterministic_model_results.Rdata'))
save(patient_diag, file = paste0(directories$dir_dat_deriv, '/deterministic_diagnostics.Rdata'))
save(outcomes_over_time_array, file = paste0(directories$dir_dat_deriv, '/deterministic_outcomes_over_time.Rdata'))
save(population, file = paste0(directories$dir_dat_deriv, '/deterministic_population.Rdata'))

################################################################################
#                                                                              #
# Cleanup                                                                      #
#                                                                              #
################################################################################

rm(model, life_tables_weight, patient_results, model_results,
   patient_diag, outcomes_over_time, outcomes_over_time_array, population)

################################################################################
#                                                                              #
# REPORT RUNTIME                                                               #
#                                                                              #
################################################################################

end.time <- Sys.time()

print(paste('Run time =', round(as.numeric(end.time, units = "secs") - as.numeric(start.time, units = "secs"), 2)/60, 'minutes', sep = ' '))

gc()