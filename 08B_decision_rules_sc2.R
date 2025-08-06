################################################################################
#                                                                              #
# Study: Cost-Effectiveness Analysis of Gene Therapy for Haemophilia B         #
# Design: Cost-Effectiveness Model using Bleeds as a continuous variable       #
# Outcome: Costs, QALYS, ICER                                                  #
# Task: Run scenario 2 for decision rules - automatic treatment success        #
# Author: Niklaus Meier                                                        #
# R version: 4.5.0                                                             #
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

nloops <- 100

## We prepare an empty data frame to save our results across multiple loops
# of our probabilistic sensitivity analysis
# We save both the Model Results Combined (MRC) and the
# Patient Results Combined (PRC)

MRC   <- data.frame()
PRC   <- data.frame()

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
  
  model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                           ncycle = (92*4)+1,
                           npatients = 100,
                           nsim = 20,
                           mode = "probabilistic",
                           treatments = c("ONDEMAND","ETRA_ONDEMAND","PROPHYLAXIS","ETRA_PROPH"),
                           threshold = 50000)
  
  model <- fun_parameter_baseline(model = model)
  model <- fun_param_sample(model = model)
  model <- fun_scenarios(model = model)
  
  ################################################################################
  #                                                                              #
  # Run Model                                                                    #
  #                                                                              #
  ################################################################################
  
  model <- fun_gen_pop(model = model, mode = "heterogeneous") # Must be heterogeneous or doesn't work right
  
  model <- fun_PROPHYLAXIS_abr(model = model)
  model <- fun_ONDEMAND_abr(model = model)
  model <- fun_ETRA_PROPH_and_OD_abr(model = model, mode = "automatic") # Very important here: Random treatment success rather than expected value makes better treatment predictions harder
  
  model <- fun_bleeds_death(model = model, mode = "expected_value") # Very important here: Random death rather than expected value makes better treatment predictions harder
  model <- fun_resources(model = model)
  model <- fun_utility(model = model)
  model <- fun_costs(model = model)
  
  ################################################################################
  #                                                                              #
  # Saving Results                                                               #
  #                                                                              #
  ################################################################################
  
  # Save results to data frame
  model_results   <- fun_results(model, mode = "simulation")
  patient_results <- fun_results(model, mode = "patient")
  
  # Combine results of this loop with results of previous loops
  MRC   <- rbind(MRC, model_results)
  PRC   <- rbind(PRC, patient_results)
  
  print(paste0(i, "/", nloops, " loops completed"))
  gc()
  
}

################################################################################
#                                                                              #
# Optimal treatment                                                            #
#                                                                              #
################################################################################

print('Identifying optimal treatment allocation')

# Identify the optimal individual treatments
dec_rules_sc2 <- fun_opt_treat(data = PRC, 
                               ntreat = length(model[["sim"]]), 
                               treatments = names(model[["sim"]]),
                               covariates = c("baseline_abr_ind", "age", "sex"),
                               optimize = "NMB")

################################################################################
#                                                                              #
# Recursive Partitioning                                                       #
#                                                                              #
################################################################################

print('Recursive Partitioning')

# Run recursive partitioning (rpart) to get optimal structure of decision tree
dec_rules_sc2 <- fun_rpart(dec_rules = dec_rules_sc2, 
                           des_split = 2,
                           treatments = names(model[["sim"]]),
                           nsim = nloops * model$struct$nsim)

################################################################################
#                                                                              #
# Save Results                                                                 #
#                                                                              #
################################################################################

save(dec_rules_sc2, file = paste0(directories$dir_dat_deriv, '/dec_rules_sc2.Rdata'))

################################################################################
#                                                                              #
# Cleanup                                                                      #
#                                                                              #
################################################################################

rm(i, life_tables_weight, nloops,
   model, model_results, MRC,
   patient_results, PRC,
   dec_rules_sc2)

################################################################################
#                                                                              #
# REPORT RUNTIME                                                               #
#                                                                              #
################################################################################

end.time <- Sys.time()

print(paste('Run time =', round(as.numeric(end.time, units = "secs") - as.numeric(start.time, units = "secs"), 2)/60, 'minutes', sep = ' '))

gc()




















