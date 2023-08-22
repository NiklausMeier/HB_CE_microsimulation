#' Function to determine utility of patients in each cycle 
#' Input: Requires a model structure that has been generated with "fun_model_setup", 
#' where ABR has been calculated based on the treatment (e.g. "fun_ETRANACOGENE_abr"),
#' and where the number and types of bleeds per cycle have been calculated with "fun_bleeds_death"

fun_utility <- function(model){ 
  
  for (i in 1:model$struct$nsim) {
    for (j in 1:length(model$sim)) {
      
      #===============================================================================
      # Baseline and Sex
      #===============================================================================
      
      # Baseline utility of 0.9508566 (Ara and Brazier 2010) is used
      # If male, utility of 0.0212126 is added
      
      model[["sim"]][[j]][[i]][,"utility_baseline_sex",] <- model[["params"]][["prob_sample"]][i,"util_baseline"] + model[["sim"]][[j]][[i]][,"sex",] * model[["params"]][["prob_sample"]][i,"util_male"]
      
      #===============================================================================
      # Age
      #===============================================================================
      
      # Per year of age, -0.0002587 disutility is applied
      # Per squared year of age, -0.0000332 additional disutility is applied
      
      model[["sim"]][[j]][[i]][,"disutility_age",] <- (model[["sim"]][[j]][[i]][,"age",] * model[["params"]][["prob_sample"]][i,"util_age"]) + ((model[["sim"]][[j]][[i]][,"age",]^2) * model[["params"]][["prob_sample"]][i,"util_agesq"])
      
      #===============================================================================
      # Arthropathy
      #===============================================================================
      
      for (k in 1:model$struct$npatients) {
        
        # We check for each patient in which cycles their Pettersson Score (PS) is between 13 and 21, or above 22
        PS_13_to_21_cycles <- which(as.vector(model[["sim"]][[j]][[i]][,"pettersson_score",k]) > 13 & as.vector(model[["sim"]][[j]][[i]][,"pettersson_score",k]) < 22 )
        PS_22_plus_cycles  <- which(as.vector(model[["sim"]][[j]][[i]][,"pettersson_score",k]) > 21)
        
        # We apply the appropriate disutility based on their PS Score
        model[["sim"]][[j]][[i]][PS_13_to_21_cycles,"disutility_arthropathy",k] <- model[["params"]][["prob_sample"]][i,"util_PS_13_to_21"]
        model[["sim"]][[j]][[i]][PS_22_plus_cycles,"disutility_arthropathy",k]  <- model[["params"]][["prob_sample"]][i,"util_PS_13_to_21"] + model[["params"]][["prob_sample"]][i,"util_PS_22_plus"]
        
      }
      
      #===============================================================================
      # Bleeds
      #===============================================================================
      
      # Apply disutility based on expected number of bleeds and type of bleed
      # Since the bleed disutility is only for a day (based on our data source) we divide by the cycle length (in days)
      
      if (model$scenarios$sc_bleed_disutility_one_day == 1) {
        disutility_joint_bleed <- (model[["sim"]][[j]][[i]][,"bleeds_joint",] * model[["params"]][["prob_sample"]][i,"util_joint_bleed"]) / model$struct$cyc2day
        disutility_other_bleed <- (model[["sim"]][[j]][[i]][,"bleeds_other",] * model[["params"]][["prob_sample"]][i,"util_other_bleed"]) / model$struct$cyc2day
        # If we assume disutility from bleed is longer, we multiply applied disutility by duration
      } else if (model$scenarios$sc_bleed_disutility_one_day == 0) {
        disutility_joint_bleed <- (model[["sim"]][[j]][[i]][,"bleeds_joint",] * model[["params"]][["prob_sample"]][i,"util_joint_bleed"] * model$scenarios$bleed_disutility_duration) / model$struct$cyc2day
        disutility_other_bleed <- (model[["sim"]][[j]][[i]][,"bleeds_other",] * model[["params"]][["prob_sample"]][i,"util_other_bleed"] * model$scenarios$bleed_disutility_duration) / model$struct$cyc2day
      }
 
      
      model[["sim"]][[j]][[i]][,"disutility_bleeds",] <- disutility_joint_bleed + disutility_other_bleed  
      
      #===============================================================================
      # Surgery
      #===============================================================================
      
      # People receive joint surgery in the first cycle their Pettersson Score (PS) reaches 28, the threshold for clinical relevance.
      # If their PS never reaches 28, they never receive surgery.
      
      for (k in 1:model$struct$npatients) {
        
        if (28 %in% as.vector(model[["sim"]][[j]][[i]][,"pettersson_score",k])) {
          
          first_PS_28_cycle <- min(which(as.vector(model[["sim"]][[j]][[i]][,"pettersson_score",k]) == model$params$pettersson$clin_rel))
          model[["sim"]][[j]][[i]][first_PS_28_cycle,"disutility_surgery",k] <- model[["params"]][["prob_sample"]][i,"util_surgery"]
          
        }
      }

      #===============================================================================
      # Factor Infusion
      #===============================================================================
      
      # We loop over the patients and calculate their prophylactic consumption of coagulation factor
      
      for (k in 1:model$struct$npatients) {
      
        # First we need to determine in which cycles patients receive prophylaxis
        # This is true for cycles in which the patients are in the treatment arms:
        #' Treatment 14: Viral Gene Therapy treatment with prophylaxis
        #' Treatment 21: Prophylaxis treatment
        
        prophylaxis_cycles <- which(as.vector(model[["sim"]][[j]][[i]][,"treatment",k]) == 14 
                                    | as.vector(model[["sim"]][[j]][[i]][,"treatment",k]) == 21)
        
        # We then apply the disutility in those cycles
        # As this is applied in terms of a "health state" utility, it is automatically scaled by the cycle length
        
        model[["sim"]][[j]][[i]][prophylaxis_cycles,"disutility_factor_infusion",k] <- model$params$prob_sample[i,"util_infusion"]
        
      }
      
      #===============================================================================
      # Total Utility
      #===============================================================================
      
      # We sum up baseline utility and all the disutilities for the total utility in this cycle
      # Total = Sex + Age + Bleeds + Arthropathy + Surgery + Factor Infusion
      
      model[["sim"]][[j]][[i]][,"utility_total",] <- model[["sim"]][[j]][[i]][,"utility_baseline_sex",] + model[["sim"]][[j]][[i]][,"disutility_age",] + model[["sim"]][[j]][[i]][,"disutility_bleeds",] + model[["sim"]][[j]][[i]][,"disutility_arthropathy",] + model[["sim"]][[j]][[i]][,"disutility_surgery",] + model[["sim"]][[j]][[i]][,"disutility_factor_infusion",]
      
      #===============================================================================
      # QALYs
      #===============================================================================
      
      # We calculate undiscounted QALYs by multiplying the undiscounted LYs by the total utility in that cycle
      
      model[["sim"]][[j]][[i]][,"QALYs_undisc",] <- model[["sim"]][[j]][[i]][,"LYs_undisc",] * model[["sim"]][[j]][[i]][,"utility_total",]
      
    }
  }
  
  return(model)
  
}