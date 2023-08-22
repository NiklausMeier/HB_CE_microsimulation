################################################################################
#                                                                              #
# Study: Cost-Effectiveness Analysis of Gene Therapy for Haemophilia B         #
# Design: Cost-Effectiveness Model using Bleeds as a continuous variable       #
# Outcome: Costs, QALYS, ICER                                                  #
# Task: Run cohort version of model with baseline parameters                   #
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

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = 1,
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
#' 2. Determine ABR based on treatment:
#' 3. Determine bleeding and death: fun_bleeds_death
#' 4. Determine resource use and costs: fun_resources
#' 5. Determine utilities and QALYs: fun_utility

#===============================================================================
# Generate Population
#===============================================================================

# We use a homogeneous population, meaning every individual has the same baseline values
model <- fun_gen_pop(model = model, mode = "homogeneous")

#===============================================================================
# Treatment
#===============================================================================

#'  In this version of the model, all treatments with ETRANACOGENE automatically succeed.
#'  The purpose of automatic success is to eliminate all randomness from the model
#'  Since it is not entirely realistic that all treatments with ETRANACOGENE would succeed,
#'  this does also change the interpretation of the model somewhat
#'  It is therefore not intended as an independent analysis, but for diagnostics and trouble-checking,
#'  and to compare a model with and without heterogeneous populations

model <- fun_ETRANACOGENE_abr(model = model, mode = "automatic")
model <- fun_PROPHYLAXIS_abr(model = model)

#===============================================================================
# Bleeding and Death
#===============================================================================

# We use the expected value of death rather than a random draw to avoid randomness.
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
# Sample Tables
#===============================================================================

table_ETRA <- data.table(model[["sim"]][["ETRANACOGENE"]][[1]][,,1])
table_PROPH <- data.table(model[["sim"]][["PROPHYLAXIS"]][[1]][,,1])

################################################################################
#                                                                              #
# Background                                                                   #
#                                                                              #
################################################################################

# We also run the model without applying treatment to receive background mortality

model <- fun_model_setup(cyc2day = settings$time$yr2day/12,
                         ncycle = 92*12,
                         npatients = 1,
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)
model <- fun_gen_pop(model = model, mode = "homogeneous")
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_utility(model = model)

table_BACK<- data.table(model[["sim"]][["ETRANACOGENE"]][[1]][,,1])

################################################################################
#                                                                              #
# Comparison of Cycle lengths                                                  #
#                                                                              #
################################################################################

# We run the entire model once with 12 month cycles, 6 month cycles, 1 month cycles, and 1 week cycles
# With this, we wish to check how the outcomes change based on cycle length

#===============================================================================
# 12 Month Cycles
#===============================================================================

model <- fun_model_setup(cyc2day = settings$time$yr2day,
                         ncycle = 92+1,
                         npatients = 1,
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)
model <- fun_gen_pop(model = model, mode = "homogeneous")
model <- fun_ETRANACOGENE_abr(model = model, mode = "automatic")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)
patient_results_cycle_12_months <- fun_results(model, mode = "simulation")

#===============================================================================
# 3 Month Cycles
#===============================================================================

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = 1,
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)
model <- fun_gen_pop(model = model, mode = "homogeneous")
model <- fun_ETRANACOGENE_abr(model = model, mode = "automatic")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)
patient_results_cycle_3_months <- fun_results(model, mode = "simulation")

#===============================================================================
# 1 Month Cycles
#===============================================================================

model <- fun_model_setup(cyc2day = settings$time$yr2day/12,
                         ncycle = (92*12)+1,
                         npatients = 1,
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)
model <- fun_gen_pop(model = model, mode = "homogeneous")
model <- fun_ETRANACOGENE_abr(model = model, mode = "automatic")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)
patient_results_cycle_1_month <- fun_results(model, mode = "simulation")

#===============================================================================
# 1 Week Cycles
#===============================================================================

model <- fun_model_setup(cyc2day = settings$time$yr2day/52,
                         ncycle = (92*52)+1,
                         npatients = 1,
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)
model <- fun_gen_pop(model = model, mode = "homogeneous")
model <- fun_ETRANACOGENE_abr(model = model, mode = "automatic")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)
patient_results_cycle_1_week <- fun_results(model, mode = "simulation")

#===============================================================================
# Comparison
#===============================================================================

# Life Years

cycle_length_comp_LYs <- data.frame(matrix(nrow = 4, ncol = 2))
colnames(cycle_length_comp_LYs) <- c("ETRANACOGENE", "PROPHYLAXIS")
rownames(cycle_length_comp_LYs) <- c("12 months", "3 months", "1 month", "1 week")

cycle_length_comp_LYs["12 months","ETRANACOGENE"] <- patient_results_cycle_12_months[,"ETRANACOGENE_LYs_disc"]
cycle_length_comp_LYs["3 months","ETRANACOGENE"]  <- patient_results_cycle_3_months[,"ETRANACOGENE_LYs_disc"]
cycle_length_comp_LYs["1 month","ETRANACOGENE"]   <- patient_results_cycle_1_month[,"ETRANACOGENE_LYs_disc"]
cycle_length_comp_LYs["1 week","ETRANACOGENE"]    <- patient_results_cycle_1_week[,"ETRANACOGENE_LYs_disc"]

cycle_length_comp_LYs["12 months","PROPHYLAXIS"]        <- patient_results_cycle_12_months[,"PROPHYLAXIS_LYs_disc"]
cycle_length_comp_LYs["3 months","PROPHYLAXIS"]         <- patient_results_cycle_3_months[,"PROPHYLAXIS_LYs_disc"]
cycle_length_comp_LYs["1 month","PROPHYLAXIS"]          <- patient_results_cycle_1_month[,"PROPHYLAXIS_LYs_disc"]
cycle_length_comp_LYs["1 week","PROPHYLAXIS"]           <- patient_results_cycle_1_week[,"PROPHYLAXIS_LYs_disc"]

cycle_length_comp_LYs

# Costs

cycle_length_comp_costs <- data.frame(matrix(nrow = 4, ncol = 2))
colnames(cycle_length_comp_costs) <- c("ETRANACOGENE", "PROPHYLAXIS")
rownames(cycle_length_comp_costs) <- c("12 months", "3 months", "1 month", "1 week")

cycle_length_comp_costs["12 months","ETRANACOGENE"] <- patient_results_cycle_12_months[,"ETRANACOGENE_cost_disc"]
cycle_length_comp_costs["3 months","ETRANACOGENE"]  <- patient_results_cycle_3_months[,"ETRANACOGENE_cost_disc"]
cycle_length_comp_costs["1 month","ETRANACOGENE"]   <- patient_results_cycle_1_month[,"ETRANACOGENE_cost_disc"]
cycle_length_comp_costs["1 week","ETRANACOGENE"]    <- patient_results_cycle_1_week[,"ETRANACOGENE_cost_disc"]

cycle_length_comp_costs["12 months","PROPHYLAXIS"]        <- patient_results_cycle_12_months[,"PROPHYLAXIS_cost_disc"]
cycle_length_comp_costs["3 months","PROPHYLAXIS"]         <- patient_results_cycle_3_months[,"PROPHYLAXIS_cost_disc"]
cycle_length_comp_costs["1 month","PROPHYLAXIS"]          <- patient_results_cycle_1_month[,"PROPHYLAXIS_cost_disc"]
cycle_length_comp_costs["1 week","PROPHYLAXIS"]           <- patient_results_cycle_1_week[,"PROPHYLAXIS_cost_disc"]

cycle_length_comp_costs

# INMB

cycle_length_INMB <- data.frame(matrix(nrow = 4, ncol = 1))
colnames(cycle_length_INMB) <- c("ETRANACOGENE vs. PROPHYLAXIS")
rownames(cycle_length_INMB) <- c("12 months", "3 months", "1 month", "1 week")

cycle_length_INMB["12 months","ETRANACOGENE vs. PROPHYLAXIS"]        <- patient_results_cycle_12_months[,"ETRANACOGENE_PROPHYLAXIS_INMB"]
cycle_length_INMB["3 months","ETRANACOGENE vs. PROPHYLAXIS"]         <- patient_results_cycle_3_months[,"ETRANACOGENE_PROPHYLAXIS_INMB"]
cycle_length_INMB["1 month","ETRANACOGENE vs. PROPHYLAXIS"]          <- patient_results_cycle_1_month[,"ETRANACOGENE_PROPHYLAXIS_INMB"]
cycle_length_INMB["1 week","ETRANACOGENE vs. PROPHYLAXIS"]           <- patient_results_cycle_1_week[,"ETRANACOGENE_PROPHYLAXIS_INMB"]

cycle_length_INMB

################################################################################
#                                                                              #
# Parameter sampling                                                           #
#                                                                              #
################################################################################

# Show example of parameter sampling:

util_infusion_distribution <- data.frame(rgamma(2000, shape = model$params$utilities$infusion$alpha, scale = model$params$utilities$infusion$beta))
colnames(util_infusion_distribution)[1] <- "utility"

ggplot(util_infusion_distribution , aes(x=utility)) + 
  geom_histogram(binwidth=0.0005, colour="black", fill="white") +
  xlab("Disutility Infusion") +
  ylab("Count") + 
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) +
  scale_x_continuous(limits=c(0,0.04),breaks=seq(0,0.04,0.01)) +
  settings$ggplot_theme

################################################################################
#                                                                              #
# Save data                                                                    #
#                                                                              #
################################################################################

save(table_ETRA, file = paste0(directories$dir_dat_deriv, '/cohort_table_ETRA.Rdata'))
save(table_PROPH, file = paste0(directories$dir_dat_deriv, '/cohort_table_PROPH.Rdata'))
save(table_BACK, file = paste0(directories$dir_dat_deriv, '/cohort_table_BACK.Rdata'))

save(cycle_length_comp_LYs, file = paste0(directories$dir_dat_deriv, '/cycle_length_comp_LYs.Rdata'))
save(cycle_length_comp_costs, file = paste0(directories$dir_dat_deriv, '/cycle_length_comp_costs.Rdata'))
save(cycle_length_INMB, file = paste0(directories$dir_dat_deriv, '/cycle_length_INMB.Rdata'))

################################################################################
#                                                                              #
# Cleanup                                                                      #
#                                                                              #
################################################################################

rm(model, life_tables_weight, 
   table_ETRA, table_PROPH, table_BACK,
   util_infusion_distribution,
   cycle_length_comp_LYs, cycle_length_comp_costs, cycle_length_INMB,
   patient_results_cycle_12_months, patient_results_cycle_3_months, patient_results_cycle_1_month, patient_results_cycle_1_week)

################################################################################
#                                                                              #
# REPORT RUNTIME                                                               #
#                                                                              #
################################################################################

end.time <- Sys.time()

print(paste('Run time =', round(as.numeric(end.time, units = "secs") - as.numeric(start.time, units = "secs"), 2)/60, 'minutes', sep = ' '))

gc()
