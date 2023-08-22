################################################################################
#                                                                              #
# Study: Cost-Effectiveness Analysis of Gene Therapy for Haemophilia B         #
# Design: Cost-Effectiveness Model using Bleeds as a continuous variable       #
# Outcome: Costs, QALYS, ICER                                                  #
# Task: Run probabilistic version of model with various parameters             #
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

library(mgcv)

#-------------------------------------------------------------------------------
# Loop
#-------------------------------------------------------------------------------

#' The probabilistic sensitivity analysis can require a lot of memory if many 
#' simulations are generated with a large number of patients.
#' Therefore, it may be necessary to run some simulations, save the outputs,
#' and then delete the model from the environment before running more simulations.
#' We do this by running multiple loops with 1 or more simulations per loop.

# Determine how many loops you want the probabilistic sensitivity analysis to run for
# Make sure that no single loop exceeds the memory of your computer

nloops <- 30

## We prepare an empty data frame to save our results across multiple loops
# of our probabilistic sensitivity analysis
model_results_combined   <- data.frame()

################################################################################
#                                                                              #
# Probabilistic Loop                                                           #
#                                                                              #
################################################################################

for (i in 1:nloops){
  
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
  
  model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                           ncycle = (92*4)+1,
                           npatients = 100,
                           nsim = 20,
                           mode = "probabilistic",
                           treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                           threshold = 50000)
  
  model <- fun_parameter_baseline(model = model)
  
  model$params$prices$ETRANACOGENE$mean   <- 2750000
  
  model <- fun_param_sample(model = model)
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
  
  model <- fun_gen_pop(model = model, mode = "heterogeneous")
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
  model_results   <- fun_results(model, mode = "simulation")
  
  # Combine results of this loop with results of previous loops
  model_results_combined   <- rbind(model_results_combined, model_results)
  
  print(paste0(i, "/", nloops, " loops completed"))
  gc()
  
}

################################################################################
#                                                                              #
# Diagnostics                                                                  #
#                                                                              #
################################################################################

# Check for convergence of results for both models and patients
# The inputs "Interval" and "Percentage" check whether result changes by more than Y% in
# an interval of X models/patients.
# IMPORTANT: When the number of simulations/patients is divided by the interval, it must lead to a whole number, not a fraction.
# GOOD: 500 patients and interval of 100. 500 / 100 = 5 (whole number)
# BAD: 500 patients and interval of 200. 500 / 200 = 2.5 (not whole number)

model_diag <- fun_diagnostics(model_results_combined, interval = 5, percentage = 0.01 )

# Convergence tables
model_diag$convergence

################################################################################
#                                                                              #
# EVPPI                                                                        #
#                                                                              #
################################################################################

inb <- model_results_combined[,"ETRANACOGENE_PROPHYLAXIS_INMB"]

EVPI <- mean(pmax(0, inb)) - max(0, mean(inb))

theta1 <- model_results_combined[,"ETRANACOGENE_max_bleed_reduction_duration"]

model_EVI <- gam(inb ~ s(theta1))
g.hat <- fitted(model_EVI)

evppi <- mean(pmax(0, g.hat)) - max(0, mean(g.hat))

plot(fitted(model_EVI), residuals(model_EVI))
gam.check(model_EVI)

# assemble prob params

param_vector  <- data.frame(cbind(model_results_combined[, "ETRANACOGENE_max_bleed_reduction_duration"],
                      model_results_combined[,  "ETRANACOGENE_abr" ]))
colnames(param_vector) <- c("ETRANACOGENE_max_bleed_reduction_duration",
                            "ETRANACOGENE_abr")

# vector of evppi variables

evppi_vector <- numeric(2)
names(evppi_vector) <- c("ETRANACOGENE_max_bleed_reduction_duration",
                         "ETRANACOGENE_abr")


for (i in (1:2)) { 
  print(i)
  parameter_of_interest <- param_vector[,i]
  model_EVI <- gam(inb ~ s(parameter_of_interest))
  fittedValues <- fitted(model_EVI)
  evppi_vector[i] <- mean(pmax(0, fittedValues)) - max(0, mean(fittedValues))
}

x.points <- barplot(
  evppi_vector, ylab = "partial EVPI", xaxt = "n", ylim = c(0, 200000)
)

axis(
  side = 1, at = x.points, labels = names(evppi_vector), 
  tick = FALSE, hadj = 0.8, las = 3
)

createInputs(param_vector)

effects <- cbind(model_results_combined[,"ETRANACOGENE_QALYs_disc"], 
                 model_results_combined[,"PROPHYLAXIS_QALYs_disc"])

costs <- cbind(model_results_combined[,"ETRANACOGENE_cost_disc"], 
               model_results_combined[,"PROPHYLAXIS_cost_disc"])

bcea.out <- bcea(eff = effects, cost = costs, ref = 1,
                 plot = T)

inputs <- createInputs(param_vector)
info.rank(bcea.out, inputs, wtp = 50000)

################################################################################
#                                                                              #
# Save Results                                                                 #
#                                                                              #
################################################################################

save(model_results_combined, file = paste0(directories$dir_dat_deriv, '/probabilistic_model_results.Rdata'))
save(model_diag, file = paste0(directories$dir_dat_deriv, '/probabilistic_model_diagnostics.Rdata'))

################################################################################
#                                                                              #
# Cleanup                                                                      #
#                                                                              #
################################################################################

rm(model, i, life_tables_weight, nloops,
   model_results, model_results_combined,
   model_diag)

################################################################################
#                                                                              #
# REPORT RUNTIME                                                               #
#                                                                              #
################################################################################

end.time <- Sys.time()

print(paste('Run time =', round(as.numeric(end.time, units = "secs") - as.numeric(start.time, units = "secs"), 2)/60, 'minutes', sep = ' '))

gc()




