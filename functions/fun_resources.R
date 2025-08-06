#' Function to determine resource use of patients in each cycle
#' We adjust for the probability of being alive and reduce resource use accordingly

fun_resources <- function(model){ 
  
  for (i in 1:model$struct$nsim) {
    for (j in 1:length(model$sim)) {
      
      #===============================================================================
      # Treatment
      #===============================================================================
      
      #-------------------------------------------------------------------------------
      # Etranacogene Dezaparvovec
      #-------------------------------------------------------------------------------
      
      # If patients are treated with viral gene therapy, treatment occurs once in first period
      # Since treatment occurs only in first period, all patients are alive
      
      if (names(model[["sim"]])[j] == "ETRANACOGENE") {
        model[["sim"]][[j]][[i]][1,"res_ETRANACOGENE",] <- 1
      }
      
      if (names(model[["sim"]])[j] == "ETRA_PROPH") {
        model[["sim"]][[j]][[i]][1,"res_ETRANACOGENE",] <- 1
      }
      
      if (names(model[["sim"]])[j] == "ETRA_ONDEMAND") {
        model[["sim"]][[j]][[i]][1,"res_ETRANACOGENE",] <- 1
      }
      
      #-------------------------------------------------------------------------------
      # Prophylaxis
      #-------------------------------------------------------------------------------
      
      # If patients are treated with prophylaxis, treatment occurs once per week
      # We first determine the number of weeks per cycle, based on cycle length
      weekspercycle <- model$struct$cyc2day/settings$time$wk2day
      
      # We loop over the patients and calculate their prophylactic consumption of coagulation factor
      
        for (k in 1:model$struct$npatients) {
          
          # First we need to determine in which cycles patients receive prophylaxis
          # This is true for cycles in which the patients are in the treatment arms:
          #' Treatment 12: Etranacogene Dezaparvovec period failure
          #' Treatment 13: Etranacogene Dezaparvovec with prophylaxis
          #' Treatment 21: Prophylaxis treatment
          
          # Check in which cycles a patient is in treatment arms 12, 14, or 21
          prophylaxis_cycles <- which(as.vector(model[["sim"]][[j]][[i]][,"treatment",k]) == 12 
                                      | as.vector(model[["sim"]][[j]][[i]][,"treatment",k]) == 14 
                                      | as.vector(model[["sim"]][[j]][[i]][,"treatment",k]) == 21) 
          
          # We determine the amount of IU needed per week for that patient by multiplying:
          # 1. Weight of patient in that cycle
          # 2. Quantity of IU/kg body weight, based on recommended dosage
          
          IU_consumption_per_week <- model[["sim"]][[j]][[i]][prophylaxis_cycles,"weight",k] * model$params$prob_sample[i,"res_IU_kg_PROPHYLAXIS"] 
          
          # We check whether vial sharing is assumed or not.
          # If not, then the IU consumption per week is rounded up to the assumed vial size.
          # We multiply:
          # 1. IU consumption per week
          # 2. Number of weeks per cycle
          # 3. Probability of being alive
          
          if (model$scenarios$sc_vial_sharing == 1) {
            
            model[["sim"]][[j]][[i]][prophylaxis_cycles,"res_IU_PROPHYLAXIS",k] <- IU_consumption_per_week * weekspercycle * model[["sim"]][[j]][[i]][prophylaxis_cycles,"alive_corrected",k]
          
            } else {
            
            model[["sim"]][[j]][[i]][prophylaxis_cycles,"res_IU_PROPHYLAXIS",k] <- round_any(IU_consumption_per_week, accuracy = model$scenarios$vial_size, f = ceiling) * weekspercycle * model[["sim"]][[j]][[i]][prophylaxis_cycles,"alive_corrected",k]
            
          }
          
        }
      
      #===============================================================================
      # Coagulation Factor to treat bleeds
      #===============================================================================
      
      # We calculate the expected quantity of needed coagulation factor to treat bleeds separately for each type of bleed
      # This is based on weight as well as the dosage for Nonacog Beta Pegol
      # We also  adjust for the probability of needing to treat an acute bleed with FIX
      
      IU_bleed_joint <- model[["sim"]][[j]][[i]][,"bleeds_joint",] * model[["sim"]][[j]][[i]][,"weight",] * model$params$prob_sample[i,"res_IU_kg_joint_bleed"] * model$params$prob_sample[i,"res_FIX_bleed_ratio"]
      IU_bleed_other <- model[["sim"]][[j]][[i]][,"bleeds_other",] * model[["sim"]][[j]][[i]][,"weight",] * model$params$prob_sample[i,"res_IU_kg_other_bleed"] * model$params$prob_sample[i,"res_FIX_bleed_ratio"]
      IU_bleed_total <- IU_bleed_joint+ IU_bleed_other
      
      # If we have perfect vial sharing, we simply use the summed up required IU of coagulation factor
      # We adjust for the probability of being alive
      
      if (model$scenarios$sc_vial_sharing == 1) {
        
        model[["sim"]][[j]][[i]][,"res_IU_treat_bleed",] <- IU_bleed_total  * model[["sim"]][[j]][[i]][,"alive_corrected",]
        
      } else {
        
        # If we do not have vial sharing, we need to adjust the vial sizes based on the expected number of bleeds
        # This is the case because we only simulate the expected number of bleeds per patient per cycle, 
        # and do not probabilistically draw the exact number of bleeds.
        # We do this by multiplying the assumed vial size with the expected number of bleeds in that cycle
        # This can lead to slightly implausible results in the individual cycles, as we assume vial sizes that do not actually exist
        # However, it assures that using expected values leads to the same result as probabilistic simulation (see "vial_sharing_adjustment.xlsx")
        
        exp_bleed_total <- model[["sim"]][[j]][[i]][,"bleeds_joint",] + model[["sim"]][[j]][[i]][,"bleeds_other",]
        
        adjusted_vial_size <- exp_bleed_total * model$scenarios$vial_size
        
        model[["sim"]][[j]][[i]][,"res_IU_treat_bleed",] <- round_any(IU_bleed_total, accuracy = adjusted_vial_size, f = ceiling)  * model[["sim"]][[j]][[i]][,"alive_corrected",]
        
      }
      
      #===============================================================================
      # Hospitalization
      #===============================================================================
      
      # To receive the expected hospital days per cycle, we multiply:
      # 1. The number of bleeds in a period
      # 2. The probability of being hospitalized per bleed
      # 3. The expected length of stay (LOS) in days per hospitalization
      # 4. The probability of being alive in that cycle
      
      hosp_days_joint <- model[["sim"]][[j]][[i]][,"bleeds_joint",] * model$params$prob_sample[i,"res_hosp_prob"] * model$params$prob_sample[i,"res_hosp_LOS"] * model[["sim"]][[j]][[i]][,"alive_corrected",]
      hosp_days_other <- model[["sim"]][[j]][[i]][,"bleeds_other",] * model$params$prob_sample[i,"res_hosp_prob"] * model$params$prob_sample[i,"res_hosp_LOS"] * model[["sim"]][[j]][[i]][,"alive_corrected",]
      
      # As all emergency bleeds are per assumption treated in the ICU, we sum these up
      model[["sim"]][[j]][[i]][,"res_hosp_ICU_days",]  <- hosp_days_joint + hosp_days_other 

      #===============================================================================
      # Surgery
      #===============================================================================
      
      # People receive joint surgery in the first cycle their Pettersson Score (PS) reaches 28, the threshold for clinical relevance.
      # If their PS never reaches 28, they never receive surgery.
      # If people receive surgery, they also spend time in the hospital ward.
      # We adjust for the probability of being alive in that cycle.
      
      for (k in 1:model$struct$npatients) {
        
        if (28 %in% as.vector(model[["sim"]][[j]][[i]][,"pettersson_score",k])) {
          first_PS_28_cycle <- min(which(as.vector(model[["sim"]][[j]][[i]][,"pettersson_score",k]) == model$params$pettersson$clin_rel))
          model[["sim"]][[j]][[i]][first_PS_28_cycle,"res_surgery",k] <- model[["sim"]][[j]][[i]][first_PS_28_cycle,"alive_corrected",k]
          
          model[["sim"]][[j]][[i]][first_PS_28_cycle,"res_hosp_ward_days",k] <- model[["sim"]][[j]][[i]][first_PS_28_cycle,"res_surgery",k] * model$params$resources$surg_LOS
        }
        
      }
      
    }
  }
      
  return(model)
  
}