################################################################################
#                                                                              #
# Study: Cost-Effectiveness Analysis of Gene Therapy for Haemophilia B         #
# Design: Cost-Effectiveness Model using Bleeds as a continuous variable       #
# Outcome: Costs, QALYS, ICER                                                  #
# Task: Use probabilistic version model as basis for individualized            #
# optimal treatment decision rules                                             #
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
  model <- fun_ETRA_PROPH_and_OD_abr(model = model, mode = "random") # Very important here: Random treatment success rather than expected value makes better treatment predictions harder
  
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
dec_rules <- fun_opt_treat(data = PRC, 
                           ntreat = length(model[["sim"]]), 
                           treatments = names(model[["sim"]]),
                           covariates = c("baseline_abr_ind", "age", "sex"),
                           optimize = "NMB")

################################################################################
#                                                                              #
# LASSO                                                                        #
#                                                                              #
################################################################################

print('LASSO')
dec_rules <- fun_LASSO(dec_rules)

################################################################################
#                                                                              #
# Recursive Partitioning                                                       #
#                                                                              #
################################################################################

print('Recursive Partitioning')

# Run recursive partitioning (rpart) to get optimal structure of decision tree
dec_rules <- fun_rpart(dec_rules = dec_rules, 
                       des_split = 2,
                       treatments = names(model[["sim"]]),
                       nsim = nloops * model$struct$nsim)

################################################################################
#                                                                              #
# POLICY TREE                                                                  #
#                                                                              #
################################################################################

print('Policy Tree')

# Apply a policy tree to identify decision rules
dec_rules <- fun_policy_tree(dec_rules = dec_rules, 
                             hybrid    = FALSE,
                             pt_depth  = 2)

# Plot policy tree
plot(dec_rules$PT$policy_tree, leaf.labels = levels(dec_rules[["W"]]))
html_plot <- plot(dec_rules$PT$policy_tree, leaf.labels = levels(dec_rules[["W"]]))
svg <- export_svg(html_plot)
rsvg::rsvg_png(charToRaw(svg), paste0(directories$dir_dat_deriv, "policy_tree.png"), width = 1300, height = 500)

################################################################################
#                                                                              #
# Convergence                                                                  #
#                                                                              #
################################################################################

# Check for convergence of results for both models and patients
# The inputs "Interval" and "Percentage" check whether result changes by more than Y% in
# an interval of X models/patients.
# IMPORTANT: When the number of simulations/patients is divided by the interval, it must lead to a whole number, not a fraction.
# GOOD: 500 patients and interval of 100. 500 / 100 = 5 (whole number)
# BAD: 500 patients and interval of 200. 500 / 200 = 2.5 (not whole number)

patient_diag <- fun_diagnostics(PRC, 
                                interval = nrow(PRC)/10, 
                                percentage = 0.01 )

# Convergence tables
patient_diag$convergence

outcomes <- c("ONDEMAND_QALYs_disc", "ETRA_ONDEMAND_QALYs_disc",
              "PROPHYLAXIS_QALYs_disc", "ETRA_PROPH_QALYs_disc",
              "ONDEMAND_cost_disc", "ETRA_ONDEMAND_cost_disc",
              "PROPHYLAXIS_cost_disc", "ETRA_PROPH_cost_disc",
              "ONDEMAND_NMB", "ETRA_ONDEMAND_NMB",
              "PROPHYLAXIS_NMB", "ETRA_PROPH_NMB")

patient_diag$plots <- fun_diagnostics_plots(patient_results[,outcomes])
patient_diag$convergence <- patient_diag$convergence[,outcomes]

## Way to take quick look at results
# dec_rules$comp_table
# count(dec_rules[["opt_treat"]][["OPTIMAL"]])
# count(dec_rules[["rpart"]][["predicted"]])

################################################################################
#                                                                              #
# Probabilistic Uncertainty with Recursive Partitioning                        #
#                                                                              #
################################################################################

# Get all possible variable combinations
dec_rules <- fun_var_combinations(dec_rules,
                                  lower_percentile = 0.05,
                                  upper_percentile = 0.95)

# Find optimal version of decision tree within each probabilistic simulation
dec_rules <- fun_rpart_prob(dec_rules)


################################################################################
#                                                                              #
# Save Results                                                                 #
#                                                                              #
################################################################################

save(MRC, file = paste0(directories$dir_dat_deriv, '/dec_rules_model_results.Rdata'))
save(PRC, file = paste0(directories$dir_dat_deriv, '/dec_rules_patient_results.Rdata'))
save(dec_rules, file = paste0(directories$dir_dat_deriv, '/dec_rules.Rdata'))
save(patient_diag, file = paste0(directories$dir_dat_deriv, '/dec_rules_diagnostics.Rdata'))

################################################################################
#                                                                              #
# Cleanup                                                                      #
#                                                                              #
################################################################################

rm(i, life_tables_weight, nloops, outcomes,
   model, model_results, MRC,
   patient_results, PRC,
   dec_rules, patient_diag)

################################################################################
#                                                                              #
# REPORT RUNTIME                                                               #
#                                                                              #
################################################################################

end.time <- Sys.time()

print(paste('Run time =', round(as.numeric(end.time, units = "secs") - as.numeric(start.time, units = "secs"), 2)/60, 'minutes', sep = ' '))

gc()
