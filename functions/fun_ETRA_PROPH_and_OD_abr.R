#' Function to determine treatment and number of bleeds per cycle after treatment with Etranacogene Dezaparvovec (ED)
#' Includes treatment switching to Prophylaxis as well as On-Demand
#' We include both treatments in the same function so that treatment success/failure is consistent in both treatment arms
#' Treatment 11: ED first treatment period success
#' Treatment 12: ED first treatment period failure
#' Treatment 13: ED observation period
#' Treatment 14: ED treatment switch to Prophylaxis
#' Treatment 31: ED treatment switch to On-Demand
#' 
#' Input:
#' 1.   Requires a model structure that has been generated with "fun_model_setup"
#' 2.   Requires a mode to be set for treatment success as "success" or "random"
#' 2.A. Under "automatic": all treatments succeed
#' 2.B. Under "random": treatment success is set to be 0 or 1 as a random chance given the probability of viral gene therapy succeeding
#' 
#'  The purpose of automatic success is to eliminate all randomness from the model
#'  Since it is not entirely realistic that all treatments with viral gene therapy would succeed,
#'  this does also change the interpretation of the model somewhat
#'  It is therefore not intended as an independent analysis, but for diagnostics and trouble-checking,
#'  and to compare a model with and without heterogeneous populations

fun_ETRA_PROPH_and_OD_abr <- function(model, mode){ 
  
  for (i in 1:model$struct$nsim) {
    
    #===============================================================================
    # Pre-determined bleed reduction vector
    #===============================================================================
    
    # We save the bleed reduction in a separate vector for later reference
    # We will apply this to the patients where treatment with viral gene therapy has succeeded later
    ETRANACOGENE_bleed_reduction_cycle = rep(0,length = model$struct$ncycle)
    
    for (j in 1:model$struct$ncycle){
      
      # Calculate how many years have passed in the current cycle
      current_cycle_years <- (j-1)*model$struct$cyc2day/settings$time$yr2day
      
      # If the current cycle has passed the point of maximum reduction, 
      # Reduce bleed reduction by percentage value specified per year past that point
      if (current_cycle_years <= model$params$prob_sample[i,"ETRANACOGENE_max_bleed_reduction_duration"]){
        
        ETRANACOGENE_bleed_reduction_cycle[j]                      <- model$params$prob_sample[i,"ETRANACOGENE_relative_bleed_reduction"]
        
      } else if (current_cycle_years > model$params$prob_sample[i,"ETRANACOGENE_max_bleed_reduction_duration"]) {
        ETRANACOGENE_bleed_reduction_cycle[j]                      <- max(model$params$prob_sample[i,"ETRANACOGENE_relative_bleed_reduction"] - ((current_cycle_years - model$params$prob_sample[i,"ETRANACOGENE_max_bleed_reduction_duration"]) * model$params$prob_sample[i,"ETRANACOGENE_bleed_increase_per_year"]),0)
        
      }
    }
    
    #===============================================================================
    # Treatment success vs failure
    #===============================================================================
    
    randomdraw <- list()
    
    if (mode == "random") {
      
      # Determine treatment success/failure:
      # Based on the probability of success, we draw a 1 (successful treatment) or a 0 (unsuccessful treatment)
      # for every patient in the ETRANACOGENE arm
      randomdraw$success <- sample(c(1,0),
                                   model$struct$npatients, 
                                   TRUE, 
                                   c(model$params$prob_sample[i,"ETRANACOGENE_success_prob"],1-model$params$prob_sample[i,"ETRANACOGENE_success_prob"]))
      
      #' Treatment 11: Etranacogene Dezaparvovec treatment period success
      #' Treatment 12: Etranacogene Dezaparvovec treatment period failure
      #' Treatment 13: Etranacogene Dezaparvovec observation period
      #' Treatment 14: Etranacogene Dezaparvovec treatment with prophylaxis
      
      ## ED INTO PROPHYLAXIS
      
      # For those patients where the treatment with ETRANACOGENE failed (treatment == 12), 
      # they receive ETRANACOGENE but also switch immediately to Prophylaxis (treatment == 14)
      model[["sim"]][["ETRA_PROPH"]][[i]][1,"treatment", randomdraw$success == 0]                       <- 12
      model[["sim"]][["ETRA_PROPH"]][[i]][2:model$struct$ncycle,"treatment", randomdraw$success == 0]   <- 14
      
      # For those patients where the treatment with ETRANACOGENE succeeded (treatment == 11),
      # they are only treated in the first cycle
      # Afterwards they are observed but not treated further, pending future treatment switches
      model[["sim"]][["ETRA_PROPH"]][[i]][1,"treatment", randomdraw$success == 1]                       <- 11
      model[["sim"]][["ETRA_PROPH"]][[i]][2:model$struct$ncycle,"treatment", randomdraw$success == 1]   <- 13
      
      ## ED INTO ON-DEMAND
      
      # For those patients where the treatment with ETRANACOGENE failed (treatment == 12), 
      # they receive ETRANACOGENE but also switch immediately to On-Demand (treatment == 31)
      model[["sim"]][["ETRA_ONDEMAND"]][[i]][1,"treatment", randomdraw$success == 0]                       <- 12
      model[["sim"]][["ETRA_ONDEMAND"]][[i]][2:model$struct$ncycle,"treatment", randomdraw$success == 0]   <- 31
      
      # For those patients where the treatment with ETRANACOGENE succeeded (treatment == 11),
      # they are only treated in the first cycle
      # Afterwards they are observed but not treated further, pending future treatment switches
      model[["sim"]][["ETRA_ONDEMAND"]][[i]][1,"treatment", randomdraw$success == 1]                       <- 11
      model[["sim"]][["ETRA_ONDEMAND"]][[i]][2:model$struct$ncycle,"treatment", randomdraw$success == 1]   <- 13
      
    } else if (mode == "automatic") {
      
      # If treatment automatically succeeds
      # Afterwards they are observed but not treated further, pending future treatment switches
      
      ## ED INTO PROPHYLAXIS
      
      model[["sim"]][["ETRA_PROPH"]][[i]][1,"treatment", ]                                              <- 11
      model[["sim"]][["ETRA_PROPH"]][[i]][2:model$struct$ncycle,"treatment", ]                          <- 13
      
      ## ED INTO ON-DEMAND
      
      model[["sim"]][["ETRA_ONDEMAND"]][[i]][1,"treatment", ]                                              <- 11
      model[["sim"]][["ETRA_ONDEMAND"]][[i]][2:model$struct$ncycle,"treatment", ]                          <- 13
      
    }
    
    failure <- (model[["sim"]][["ETRA_PROPH"]][[i]][1,"treatment",] == 12)
    
    #===============================================================================
    # Apply bleed reduction and calculate ABR
    #===============================================================================
    
    ## ED INTO PROPHYLAXIS
    
    # For all patients where treatment failed, we know they will be receiving Prophylaxis for the rest of the model
    
    model[["sim"]][["ETRA_PROPH"]][[i]][,"abr_reduction", failure == 1] <- model$params$prob_sample[i,"PROPHYLAXIS_relative_bleed_reduction"]
    model[["sim"]][["ETRA_PROPH"]][[i]][,"current_abr", failure == 1]   <- model[["sim"]][["ETRA_PROPH"]][[i]][,"baseline_abr",failure == 1] * (1 - model[["sim"]][["ETRA_PROPH"]][[i]][,"abr_reduction",failure == 1])
    
    # For patients where treatment succeeded
    
    model[["sim"]][["ETRA_PROPH"]][[i]][,"abr_reduction",failure == 0]  <- ETRANACOGENE_bleed_reduction_cycle
    model[["sim"]][["ETRA_PROPH"]][[i]][,"current_abr",failure == 0]    <- model[["sim"]][["ETRA_PROPH"]][[i]][,"baseline_abr",failure == 0] * (1 - model[["sim"]][["ETRA_PROPH"]][[i]][,"abr_reduction",failure == 0])
    
    ## ED INTO ON-DEMAND
    
    # For all patients where treatment failed, we know they will be receiving On-Demand for the rest of the model
    
    model[["sim"]][["ETRA_ONDEMAND"]][[i]][,"abr_reduction", failure == 1] <- 0
    model[["sim"]][["ETRA_ONDEMAND"]][[i]][,"current_abr",failure == 1]    <- model[["sim"]][["ETRA_ONDEMAND"]][[i]][,"baseline_abr",failure == 1] * (1 - model[["sim"]][["ETRA_ONDEMAND"]][[i]][,"abr_reduction",failure == 1])
    
    # For patients where treatment succeeded
    
    model[["sim"]][["ETRA_ONDEMAND"]][[i]][,"abr_reduction",failure == 0] <- ETRANACOGENE_bleed_reduction_cycle
    model[["sim"]][["ETRA_ONDEMAND"]][[i]][,"current_abr",failure == 0]   <- model[["sim"]][["ETRA_ONDEMAND"]][[i]][,"baseline_abr",failure == 0] * (1 - model[["sim"]][["ETRA_ONDEMAND"]][[i]][,"abr_reduction",failure == 0])
    
    #===============================================================================
    # Treatment Switching
    #===============================================================================
    
    ## ED INTO PROPHYLAXIS
    
    if (model[["scenarios"]][["sc_treatment_switching"]] == 1) {
      
      for (j in 1:model$struct$npatients) {
        
        # We check if treatment succeeded initially (treatment == 11) and whether treatment switching is even relevant
        if (model[["sim"]][["ETRA_PROPH"]][[i]][1,"treatment", j] == 11) {
          
          # Check ETRA_PROPH and treatment switch for cycles where current ABR exceeds treatment_switch_abr
          treatment_switch_cycles <- which(as.vector(model[["sim"]][["ETRA_PROPH"]][[i]][,"current_abr",j]) > model$scenarios$treatment_switch_abr)
          model[["sim"]][["ETRA_PROPH"]][[i]][treatment_switch_cycles,"treatment", j]  <- 14
          
          # Determine current bleed reduction for Prophylaxis due to treatment switching (treatment == 14):
          # Formula: 1 - ((1 - Prophylaxis Bleed Reduction)*(1 - ETRANACOGENE bleed reduction))
          # This creates a vector the length of which is all cycles in which the patient switched treatment
          # This vector overwrites the previous bleed reduction from ETRANACOGENE alone
          model[["sim"]][["ETRA_PROPH"]][[i]][treatment_switch_cycles,"abr_reduction",j] <- 1 - (1 - model$params$prob_sample[i,"PROPHYLAXIS_relative_bleed_reduction"]) * (1 - ETRANACOGENE_bleed_reduction_cycle[treatment_switch_cycles])
          model[["sim"]][["ETRA_PROPH"]][[i]][,"current_abr",j]   <- model[["sim"]][["ETRA_PROPH"]][[i]][,"baseline_abr",j] * (1 - model[["sim"]][["ETRA_PROPH"]][[i]][,"abr_reduction",j])
          
        }
      }
    }
    
    ## ED INTO ON-DEMAND
    # Patients always receive on-demand, so no need for treatment switching
    
  }
  
  return(model)
  
}
