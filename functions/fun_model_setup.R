#' Function to setup model structure
#' Inputs:
#' cyc2day: Number of days per day
#' ncycle: Number of cycles per simulation
#' npatients: Number of patients per simulation
#' nsim: Number of simulations
#' treatments: Vector of treatment names

fun_model_setup <- function(cyc2day, ncycle, npatients, nsim, mode, treatments, threshold){ 
  
  ################################################################################
  #                                                                              #
  #  Model Structure                                                             #
  #                                                                              #
  ################################################################################
  
  model <- list()
  
  # Choose cycle length, number of cycles, discount rate and the number of patients
  # Cycle length should be no shorter than 27 days, due to the calculation of disutilities for some events
  
  model$struct$cyc2day          <- cyc2day
  model$struct$ncycle           <- ncycle
  model$struct$npatients        <- npatients
  model$struct$nsim             <- nsim
  model$struct$mode             <- mode
  model$struct$threshold        <- threshold
  
  ################################################################################
  #                                                                              #
  #  Measures                                                                    #
  #                                                                              #
  ################################################################################
  
  ## Set up measures that will be used, tracked and calculated in the model array
  
  model$measures    <- c("sex", "age", "weight", "alive", "alive_corrected" , "treatment",
                         
                         "baseline_abr", "abr_reduction", "current_abr",
                         "bleeds_joint", "bleeds_other", 
                         "bleeds_joint_cum", "pettersson_score",
                         
                         "bleed_mortality_ratio", "death_prob_age", "death_prob_total",
                         
                         "utility_baseline_sex", "disutility_age", "disutility_bleeds",
                         "disutility_arthropathy", "disutility_surgery", "disutility_factor_infusion",
                         "utility_total",
                         
                         "res_ETRANACOGENE", "res_IU_PROPHYLAXIS",
                         "res_IU_treat_bleed", "res_hosp_ward_days", "res_hosp_ICU_days", 
                         "res_surgery",
                         
                         "cost_ETRANACOGENE_undisc", "cost_PROPHYLAXIS_undisc", 
                         "cost_factor_treat_bleed_undisc", "cost_hospitalization_undisc", "cost_surgery_undisc",
                         
                         "cost_undisc", "LYs_undisc", "QALYs_undisc")
  
  ################################################################################
  #                                                                              #
  #  Generate simulation list with n patients per simulation                     #
  #                                                                              #
  ################################################################################
  
  # We create a 3-dimensional array for each simulation, one per treatment arm.
  # Dimension 1(rows): Cycles of time
  # Dimension 2(columns): Measures that we are tracking for the model for every patient
  # Dimension 3(depth): Patients
  # Example: Model with 100 cycles, 26 measures, and 500 patients will be a 3-dimensional array with the dimensions 100 x 26 x 500
  # For the probabilistic model, each simulation gets its own 3-dimensional array
  
  model[["sim"]] <- vector(mode = "list", length = length(treatments))
  
  names(model[["sim"]]) <- treatments
  
  # We create a 3-dimensional array for each simulation, one per treatment arm.
  # Dimension 1(rows): Cycles of time
  # Dimension 2(columns): Measures that we are tracking for the model for every patient
  # Dimension 3(depth): Patients
  # Example: Model with 100 cycles, 26 measures, and 500 patients will be a 3-dimensional array with the dimensions 100 x 26 x 500
  # For the probabilistic model, each simulation gets its own 3-dimensional array
  
  for (i in 1:model$struct$nsim) {
    for (j in 1:length(treatments)) {
      
      model[["sim"]][[j]][[i]]         <- array(0,
                                                dim=c(model$struct$ncycle,length(model$measures),model$struct$npatients),
                                                dimnames = list(1:model$struct$ncycle,
                                                                model$measures,
                                                                1:model$struct$npatients))
    }
  }
  
  return(model)
  
}