################################################################################
#                                                                              #
# Study: Cost-Effectiveness Analysis of Gene Therapy for Haemophilia B         #
# Design: Cost-Effectiveness Model using Bleeds as a continuous variable       #
# Outcome: Costs, QALYS, ICER                                                  #
# Task: Analyze scenarios                                                      #
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

# We load the population from our deterministic analysis to use it here as well
load(file = paste0(directories$dir_dat_deriv, "/deterministic_population.Rdata"))

################################################################################
#                                                                              #
# Scenario Parameters and Ranges                                               #
#                                                                              #
################################################################################

scen_analysis <- list()

#===============================================================================
# Parameters
#===============================================================================

# We select and name the parameters for which we want to perform our scenarios
# These names NEED to correspond exactly to those in fun_param_sample

#-------------------------------------------------------------------------------
# ETRANACOGENE
#-------------------------------------------------------------------------------

scen_analysis$ETRANACOGENE_relative_bleed_reduction$name           <- "ETRANACOGENE_relative_bleed_reduction"
scen_analysis$ETRANACOGENE_max_bleed_reduction_duration$name       <- "ETRANACOGENE_max_bleed_reduction_duration"
scen_analysis$ETRANACOGENE_bleed_increase_per_year$name            <- "ETRANACOGENE_bleed_increase_per_year"
scen_analysis$ETRANACOGENE_success_prob$name                       <- "ETRANACOGENE_success_prob"

#-------------------------------------------------------------------------------
# Prices
#-------------------------------------------------------------------------------

scen_analysis$price_coagulation_factor_IU$name                  <- "price_coagulation_factor_IU"
scen_analysis$price_ETRANACOGENE$name                           <- "price_ETRANACOGENE"
scen_analysis$price_hosp_ICU_days$name                          <- "price_hosp_ICU_days"

#===============================================================================
# Value Ranges
#===============================================================================

# We define different ranges of values for our parameters, to see how costs and cost-effectiveness changes based on this

#-------------------------------------------------------------------------------
# ETRANACOGENE
#-------------------------------------------------------------------------------

scen_analysis$ETRANACOGENE_relative_bleed_reduction$range           <- seq(0.1, 1.0, by = 0.1)
scen_analysis$ETRANACOGENE_max_bleed_reduction_duration$range       <- c(1,5,10,15,20,25,30)
scen_analysis$ETRANACOGENE_bleed_increase_per_year$range            <- seq(0.05, 0.50, by = 0.05)
scen_analysis$ETRANACOGENE_success_prob$range                       <- seq(0.5, 1.0)

#-------------------------------------------------------------------------------
# Prices
#-------------------------------------------------------------------------------

scen_analysis$price_coagulation_factor_IU$range                  <- seq(0, 3, by = 0.5)
scen_analysis$price_ETRANACOGENE$range                           <- seq(1000000, 3000000, by = 500000)
scen_analysis$price_hosp_ICU_days$range                          <- 6000

################################################################################
#                                                                              #
# Scenario Analysis #1                                                         #
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Prepare Model
#-------------------------------------------------------------------------------

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = length(population$ETRANACOGENE[[1]][1,1,]),
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE","PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)

# We save our baseline parameters
param_sample <- model$params$prob_sample

#-------------------------------------------------------------------------------
# Loop
#-------------------------------------------------------------------------------

#' Testing scenarios can require a lot of memory if many 
#' simulations are generated with a large number of patients.
#' Therefore, it may be necessary to run some simulations, save the outputs,
#' and then delete the model from the environment before running more simulations.
#' We do this by running multiple loops with 1 simulation per loop.

## We prepare an empty data frame to save our results across multiple loops

model_results_combined   <- data.frame()
patient_results_combined <- data.frame()

for (j in 1:length(scen_analysis)) {
  
  print(paste0("Parameter ", j, "/", length(scen_analysis), ": " ,  scen_analysis[[j]][["name"]]))
  
  for (i in 1:length(scen_analysis[[j]][["range"]])) {
    
    # We reset the population in our model to that of our baseline model
    
    model$sim <- population
    
    # We reset parameters in model
    
    model$params$prob_sample <- param_sample
    
    # We assign the correct parameter value for our scenario analysis
    
    model[["params"]][["prob_sample"]][,paste0(scen_analysis[[j]][["name"]])] <- scen_analysis[[j]][["range"]][[i]]
    
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
    
    patient_results[,"scen_analysis_param"] <- scen_analysis[[j]][["name"]]
    model_results[,"scen_analysis_param"]  <- scen_analysis[[j]][["name"]]
    
    # Combine results of this loop with results of previous loops
    patient_results_combined <- rbind(patient_results_combined, patient_results)
    model_results_combined   <- rbind(model_results_combined, model_results)
    
    # Print progress
    print(paste0(i, "/", length(scen_analysis[[j]][["range"]]), " simulations complete"))
    
    # Clean up memory
    gc()
  }
  
}

# Save scenario analysis results to list based on parameter

for (j in 1:length(scen_analysis)) {
  scen_analysis[[j]][["patient_results"]] <- subset(patient_results_combined,scen_analysis_param == paste0(scen_analysis[[j]][["name"]]))
  scen_analysis[[j]][["model_results"]]   <- subset(model_results_combined,scen_analysis_param == paste0(scen_analysis[[j]][["name"]]))
}

################################################################################
#                                                                              #
# Scenario Analysis #2                                                         #
#                                                                              #
################################################################################

#-------------------------------------------------------------------------------
# Base case
#-------------------------------------------------------------------------------

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = length( population$ETRANACOGENE[[1]][1,1,]),
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)

# We reset the population in our model to that of our baseline model
model$sim <- population

model <- fun_ETRANACOGENE_abr(model = model, mode = "random")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)

# Save results of model to data frame
patient_results <- fun_results(model, mode = "patient")
model_results   <- fun_results(model, mode = "simulation")

patient_results[,"scen_analysis_param"] <- "base_case"
model_results[,"scen_analysis_param"]   <- "base_case"

# Combine results of this loop with results of previous loops
patient_results_combined <- rbind(patient_results_combined, patient_results)
model_results_combined   <- rbind(model_results_combined, model_results)

# Add results to scenario analysis list
scen_analysis$base_case$name <- "base_case"
scen_analysis$base_case$patient_results <- subset(patient_results_combined, scen_analysis_param == "base_case")
scen_analysis$base_case$model_results   <- subset(model_results_combined, scen_analysis_param == "base_case")

# We also test some discrete scenarios based on altering fun_scenarios

#-------------------------------------------------------------------------------
# Longer disutility from bleed
#-------------------------------------------------------------------------------

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = length(population$ETRANACOGENE[[1]][1,1,]),
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)

# We assume bleed lasts longer than one day, instead a week
model$scenarios$sc_bleed_disutility_one_day <- 0
model$scenarios$bleed_disutility_duration   <- 7

# We reset the population in our model to that of our baseline model
model$sim <- population

model <- fun_ETRANACOGENE_abr(model = model, mode = "random")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)

# Save results of model to data frame
patient_results <- fun_results(model, mode = "patient")
model_results   <- fun_results(model, mode = "simulation")

patient_results[,"scen_analysis_param"] <- "scen_bleed_disutility"
model_results[,"scen_analysis_param"]  <- "scen_bleed_disutility"

# Combine results of this loop with results of previous loops
patient_results_combined <- rbind(patient_results_combined, patient_results)
model_results_combined   <- rbind(model_results_combined, model_results)

# Add results to scenario analysis list
scen_analysis$scen_bleed_disutility$name <- "scen_bleed_disutility"
scen_analysis$scen_bleed_disutility$patient_results <- subset(patient_results_combined, scen_analysis_param == "scen_bleed_disutility")
scen_analysis$scen_bleed_disutility$model_results   <- subset(model_results_combined, scen_analysis_param == "scen_bleed_disutility")

#-------------------------------------------------------------------------------
# No treatment switching for Etranacogene Dezaparvovec
#-------------------------------------------------------------------------------

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = length( population$ETRANACOGENE[[1]][1,1,]),
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)

# We assume bleed lasts longer than one day, instead a week
model$scenarios$sc_treatment_switching <- 0

# We reset the population in our model to that of our baseline model
model$sim <- population

model <- fun_ETRANACOGENE_abr(model = model, mode = "random")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)

# Save results of model to data frame
patient_results <- fun_results(model, mode = "patient")
model_results   <- fun_results(model, mode = "simulation")

patient_results[,"scen_analysis_param"] <- "scen_ETRANACOGENE_no_switch"
model_results[,"scen_analysis_param"]   <- "scen_ETRANACOGENE_no_switch"

# Combine results of this loop with results of previous loops
patient_results_combined <- rbind(patient_results_combined, patient_results)
model_results_combined   <- rbind(model_results_combined, model_results)

# Add results to scenario analysis list
scen_analysis$scen_ETRANACOGENE_no_switch$name <- "scen_ETRANACOGENE_no_switch"
scen_analysis$scen_ETRANACOGENE_no_switch$patient_results <- subset(patient_results_combined, scen_analysis_param == "scen_ETRANACOGENE_no_switch")
scen_analysis$scen_ETRANACOGENE_no_switch$model_results   <- subset(model_results_combined, scen_analysis_param == "scen_ETRANACOGENE_no_switch")

#-------------------------------------------------------------------------------
# No vial sharing FIX
#-------------------------------------------------------------------------------

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = length( population$ETRANACOGENE[[1]][1,1,]),
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)

# We assume bleed lasts longer than one day, instead a week
model$scenarios$sc_vial_sharing <- 0

# We reset the population in our model to that of our baseline model
model$sim <- population

model <- fun_ETRANACOGENE_abr(model = model, mode = "random")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)

# Save results of model to data frame
patient_results <- fun_results(model, mode = "patient")
model_results   <- fun_results(model, mode = "simulation")

patient_results[,"scen_analysis_param"] <- "scen_no_vial_sharing"
model_results[,"scen_analysis_param"]   <- "scen_no_vial_sharing"

# Combine results of this loop with results of previous loops
patient_results_combined <- rbind(patient_results_combined, patient_results)
model_results_combined   <- rbind(model_results_combined, model_results)

# Add results to scenario analysis list
scen_analysis$scen_no_vial_sharing$name <- "scen_no_vial_sharing"
scen_analysis$scen_no_vial_sharing$patient_results <- subset(patient_results_combined, scen_analysis_param == "scen_no_vial_sharing")
scen_analysis$scen_no_vial_sharing$model_results   <- subset(model_results_combined, scen_analysis_param == "scen_no_vial_sharing")

#-------------------------------------------------------------------------------
# Starting age of 18
#-------------------------------------------------------------------------------

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = length( population$ETRANACOGENE[[1]][1,1,]),
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)

# We reset the population in our model to that of our baseline model
model$sim <- population

# Recalculate age of everyone in model to start at age 18
# First we create a vector of ages suitable for our time horizon and cycle length

age_18 <- vector(length=model$struct$ncycle)

for (i in 1:model$struct$ncycle) {
  age_18[i] <- 18 + (i-1) * model$struct$cyc2day/settings$time$yr2day
}

# We manually overwrite the age vector of each patient with our new vector

for (i in 1:length(model$sim)) {
  for (j in 1:model$struct$npatients) {
    
    model[["sim"]][[i]][[1]][,"age",j] <- age_18
    
  }
}

# We recalculate weight
# Create weight array

weight_18 <- array(0,dim = c(model$struct$ncycle,model$struct$npatients))

# Determine weight based on age and sex

for (j in 1:model$struct$npatients) {
  if (model[["sim"]][[i]][[1]][1,"sex",j] == 1) {
    weight_18[,j] <- lookup(pmin(floor(age_18),110), 
                              model$life_tables_weight$age, 
                              model$life_tables_weight$male_weight_kg)
  } else if (model[["sim"]][[i]][[1]][1,"sex",j] == 0) {
    weight_18[,j] <- lookup(pmin(floor(age_18),110), 
                              model$life_tables_weight$age, 
                              model$life_tables_weight$female_weight_kg)
  }
}

# Assign weight to patients by age and sex in each cycle
for (j in 1:length(model$sim)) {
  model[["sim"]][[j]][[1]][,"weight",]          <- weight_18
}

# We assign background mortality based on Age
# Create death-rate-age Array

death_rate_age_18 <- array(0,dim = c(model$struct$ncycle,model$struct$npatients))

# Determine rate of death based on age and sex
# We do this inside the population-generation function since it needs to be appropriate to the age of each patient in the model.

for (j in 1:model$struct$npatients) {
  if (model[["sim"]][[i]][[1]][1,"sex",j] == 1) {
    death_rate_age_18[,j] <- lookup(pmin(floor(age_18), 110), 
                                            model$life_tables_weight$age, 
                                            model$life_tables_weight$male_yearly_death_rate)
  } else if (model[["sim"]][[i]][[1]][1,"sex",j] == 0) {
    death_rate_age_18[,j] <- lookup(pmin(floor(age_18), 110), 
                                            model$life_tables_weight$age, 
                                            model$life_tables_weight$female_yearly_death_rate)
  }
}

# Create death-probability-age Array

death_prob_age_18 <- array(0,dim = c(model$struct$ncycle,model$struct$npatients))

# Convert rate into probability, adjusted for model cycle length

death_prob_age_18 <- (1 - exp(-death_rate_age_18 * (model$struct$cyc2day/settings$time$yr2day)))

# Assign probability of death based on age to patients in each cycle

for (j in 1:length(model$sim)) {
  model[["sim"]][[j]][[1]][,"death_prob_age",]          <- death_prob_age_18
}

# We also need to recalculate the number of cumulative joint bleeds in the first cycle 

for (i in 1:length(model$sim)) {
  for (j in 1:model$struct$npatients) {
    
    # Recalculate historic ABR of patient
    hist_abr <- model[["sim"]][[i]][[1]][1,"baseline_abr",j] * model$params$prob_sample[1,"hist_baseline_abr_ratio"]
    
    model[["sim"]][[i]][[1]][1,"bleeds_joint_cum",j] <- model[["sim"]][[i]][[1]][1,"age",j] * hist_abr * model$params$prob_sample[1,"share_joint_bleeds"]
    
  }
}

model <- fun_ETRANACOGENE_abr(model = model, mode = "random")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)

# Save results of model to data frame
patient_results <- fun_results(model, mode = "patient")
model_results   <- fun_results(model, mode = "simulation")

patient_results[,"scen_analysis_param"] <- "scen_age_18"
model_results[,"scen_analysis_param"]   <- "scen_age_18"

# Combine results of this loop with results of previous loops
patient_results_combined <- rbind(patient_results_combined, patient_results)
model_results_combined   <- rbind(model_results_combined, model_results)

# Add results to scenario analysis list
scen_analysis$age_18$name <- "scen_age_18"
scen_analysis$age_18$patient_results <- subset(patient_results_combined, scen_analysis_param == "scen_age_18")
scen_analysis$age_18$model_results   <- subset(model_results_combined, scen_analysis_param == "scen_age_18")

#-------------------------------------------------------------------------------
# Starting age of 60
#-------------------------------------------------------------------------------

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = length( population$ETRANACOGENE[[1]][1,1,]),
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)

# We reset the population in our model to that of our baseline model
model$sim <- population

# Recalculate age of everyone in model to start at age 60
# First we create a vector of ages suitable for our time horizon and cycle length

age_60 <- vector(length=model$struct$ncycle)

for (i in 1:model$struct$ncycle) {
  age_60[i] <- 60 + (i-1) * model$struct$cyc2day/settings$time$yr2day
}

# We manually overwrite the age vector of each patient with our new vector

for (i in 1:length(model$sim)) {
  for (j in 1:model$struct$npatients) {
    
    model[["sim"]][[i]][[1]][,"age",j] <- age_60
    
  }
}


# We recalculate weight
# Create weight array

weight_60 <- array(0,dim = c(model$struct$ncycle,model$struct$npatients))

# Determine weight based on age and sex

for (j in 1:model$struct$npatients) {
  if (model[["sim"]][[i]][[1]][1,"sex",j] == 1) {
    weight_60[,j] <- lookup(pmin(floor(age_60),110), 
                            model$life_tables_weight$age, 
                            model$life_tables_weight$male_weight_kg)
  } else if (model[["sim"]][[i]][[1]][1,"sex",j] == 0) {
    weight_60[,j] <- lookup(pmin(floor(age_60),110), 
                            model$life_tables_weight$age, 
                            model$life_tables_weight$female_weight_kg)
  }
}

# Assign weight to patients by age and sex in each cycle
for (j in 1:length(model$sim)) {
  model[["sim"]][[j]][[1]][,"weight",]          <- weight_60
}

# We assign background mortality based on Age
# Create death-rate-age Array

death_rate_age_60 <- array(0,dim = c(model$struct$ncycle,model$struct$npatients))

# Determine rate of death based on age and sex
# We do this inside the population-generation function since it needs to be appropriate to the age of each patient in the model.

for (j in 1:model$struct$npatients) {
  if (model[["sim"]][[i]][[1]][1,"sex",j] == 1) {
    death_rate_age_60[,j] <- lookup(pmin(floor(age_60), 110), 
                                    model$life_tables_weight$age, 
                                    model$life_tables_weight$male_yearly_death_rate)
  } else if (model[["sim"]][[i]][[1]][1,"sex",j] == 0) {
    death_rate_age_60[,j] <- lookup(pmin(floor(age_60), 110), 
                                    model$life_tables_weight$age, 
                                    model$life_tables_weight$female_yearly_death_rate)
  }
}

# Create death-probability-age Array

death_prob_age_60 <- array(0,dim = c(model$struct$ncycle,model$struct$npatients))

# Convert rate into probability, adjusted for model cycle length

death_prob_age_60 <- (1 - exp(-death_rate_age_60 * (model$struct$cyc2day/settings$time$yr2day)))

# Assign probability of death based on age to patients in each cycle

for (j in 1:length(model$sim)) {
  model[["sim"]][[j]][[1]][,"death_prob_age",]          <- death_prob_age_60
}


# We also need to recalculate the number of cumulative joint bleeds in the first cycle 

for (i in 1:length(model$sim)) {
  for (j in 1:model$struct$npatients) {
    
    # Recalculate historic ABR of patient
    hist_abr <- model[["sim"]][[i]][[1]][1,"baseline_abr",j] * model$params$prob_sample[1,"hist_baseline_abr_ratio"]
    
    model[["sim"]][[i]][[1]][1,"bleeds_joint_cum",j] <- model[["sim"]][[i]][[1]][1,"age",j] * hist_abr * model$params$prob_sample[1,"share_joint_bleeds"]
    
  }
}

model <- fun_ETRANACOGENE_abr(model = model, mode = "random")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)
model <- fun_utility(model = model)
model <- fun_costs(model = model)

# Save results of model to data frame
patient_results <- fun_results(model, mode = "patient")
model_results   <- fun_results(model, mode = "simulation")

patient_results[,"scen_analysis_param"] <- "scen_age_60"
model_results[,"scen_analysis_param"]   <- "scen_age_60"

# Combine results of this loop with results of previous loops
patient_results_combined <- rbind(patient_results_combined, patient_results)
model_results_combined   <- rbind(model_results_combined, model_results)

# Add results to scenario analysis list
scen_analysis$age_60$name <- "scen_age_60"
scen_analysis$age_60$patient_results <- subset(patient_results_combined, scen_analysis_param == "scen_age_60")
scen_analysis$age_60$model_results   <- subset(model_results_combined, scen_analysis_param == "scen_age_60")

#-------------------------------------------------------------------------------
# German Utilities from Grochtendreis et al. (2019)
#-------------------------------------------------------------------------------

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,
                         ncycle = (92*4)+1,
                         npatients = length(population$ETRANACOGENE[[1]][1,1,]),
                         nsim = 1,
                         mode = "deterministic",
                         treatments = c("ETRANACOGENE", "PROPHYLAXIS"),
                         threshold = 50000)

model <- fun_parameter_baseline(model = model)
model <- fun_param_sample(model = model)
model <- fun_scenarios(model = model)

# We reset the population in our model to that of our baseline model
model$sim <- population

model <- fun_ETRANACOGENE_abr(model = model, mode = "random")
model <- fun_PROPHYLAXIS_abr(model = model)
model <- fun_bleeds_death(model = model, mode = "expected_value")
model <- fun_resources(model = model)

# We use different utility function for this scenario

model <- fun_utility_scen_Grochtendreis(model = model)
model <- fun_costs(model = model)

# Save results of model to data frame
patient_results <- fun_results(model, mode = "patient")
model_results   <- fun_results(model, mode = "simulation")

patient_results[,"scen_analysis_param"] <- "scen_util_german"
model_results[,"scen_analysis_param"]  <- "scen_util_german"

# Combine results of this loop with results of previous loops
patient_results_combined <- rbind(patient_results_combined, patient_results)
model_results_combined   <- rbind(model_results_combined, model_results)

# Add results to scenario analysis list
scen_analysis$scen_util_german$name <- "scen_util_german"
scen_analysis$scen_util_german$patient_results <- subset(patient_results_combined, scen_analysis_param == "scen_util_german")
scen_analysis$scen_util_german$model_results   <- subset(model_results_combined, scen_analysis_param == "scen_util_german")

################################################################################
#                                                                              #
# Save Results                                                                 #
#                                                                              #
################################################################################

save(patient_results_combined, file = paste0(directories$dir_dat_deriv, '/scenarios_patients.Rdata'))
save(model_results_combined, file = paste0(directories$dir_dat_deriv, '/scenarios_model.Rdata'))
save(scen_analysis, file = paste0(directories$dir_dat_deriv, '/scenarios_inputs_results.Rdata'))

################################################################################
#                                                                              #
# Cleanup                                                                      #
#                                                                              #
################################################################################

rm(model, i, life_tables_weight,
   model_results, patient_results,
   param_sample, population,
   hist_abr, age_18, age_60, weight_18, weight_60, 
   death_rate_age_18, death_rate_age_60, death_prob_age_18, death_prob_age_60,
   patient_results_combined, model_results_combined, scen_analysis)

################################################################################
#                                                                              #
# REPORT RUNTIME                                                               #
#                                                                              #
################################################################################

end.time <- Sys.time()

print(paste('Run time =', round(as.numeric(end.time, units = "secs") - as.numeric(start.time, units = "secs"), 2)/60, 'minutes', sep = ' '))

gc()

