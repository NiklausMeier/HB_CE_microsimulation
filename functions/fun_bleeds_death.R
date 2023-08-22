#' Function to determine how many bleeds of each kind are experienced by a patient each cycle
#' Based on the number of joint bleeds, these are accumulated and lead to an increased Pettersson Score
#' Based on the number of total bleeds, age-standard mortality is adjusted upward
#' Input: 
#' 1. Requires a model structure that has been generated with "fun_model_setup", and where ABR has been calculated based on the treatment (e.g. "fun_ETRANACOGENE_abr")
#' 2. Requires a mode to be set for mortality as "expected_value" or "random"
#' 2.A. Under "expected_value": the alive variable takes the probability of a patient being alive at the start of a cycle
#' 2.B. Under "random": alive is set to be 0 or 1 as a random chance given the probability of living from last cycle to the current one


fun_bleeds_death <- function(model, mode){ 
  
  for (i in 1:model$struct$nsim) {
    for (j in 1:length(model$sim)) {
      
      #===============================================================================
      # Bleeds per cycle by type
      #===============================================================================
      
      for (k in 1:model$struct$npatients) {
        
        # We multiply the current ABR by the share of bleeds in each bleed category, and by the cycle length (as fraction of a year)
        # The share of "other" bleeds is simply 1 minus the other categories
        
        model[["sim"]][[j]][[i]][,"bleeds_joint",k] <-  model[["sim"]][[j]][[i]][,"current_abr",k] * model[["params"]][["prob_sample"]][i,"share_joint_bleeds"] * (model$struct$cyc2day/settings$time$yr2day)
        model[["sim"]][[j]][[i]][,"bleeds_other",k] <-  model[["sim"]][[j]][[i]][,"current_abr",k] * model[["params"]][["prob_sample"]][i,"share_other_bleeds"] * (model$struct$cyc2day/settings$time$yr2day)
        
      }
      
      #===============================================================================
      # Cumulative joint bleeds
      #===============================================================================
      
      # Joint bleeds accumulate over time, in all treatment arms
      # As this accumulates from cycle to cycle, we loop over the cycles (starting with cycle 2) chronologically, rather than over the patients
      
      for (k in 2:model$struct$ncycle) {
        
        model[["sim"]][[j]][[i]][k,"bleeds_joint_cum",] <- model[["sim"]][[j]][[i]][k-1,"bleeds_joint",] + model[["sim"]][[j]][[i]][k-1,"bleeds_joint_cum",]
        
      }
      
      #===============================================================================
      # Pettersson Score
      #===============================================================================
      
      # Based on current number of accumulated bleeds (generally accumulated bleeds divided by 12.6)
      # A Pettersson score for arthropathy is determined
      # This is later used to determine disutility from arthropathy
      
      model[["sim"]][[j]][[i]][,"pettersson_score",] <- pmin(floor(model[["sim"]][[j]][[i]][,"bleeds_joint_cum",] / model[["params"]][["prob_sample"]][i,"PS_bleeds"]), 
                                                             model$params$pettersson$max)
      
      #===============================================================================
      # Standardized Mortality Ratio (SMR) adjusted for bleeding
      #===============================================================================
      
      # We compare current ABR with the historic ABR in each cycle for all patients.
      # We calculate what the ratio of current ABR is relative to historic ABR (temp_bleed_ratio).
      # Based on this reduction, we also reduce the Standardized Mortality Ratio (SMR) of hemophilia by the same percentage.
      
      temp_bleed_ratio <- model[["sim"]][[j]][[i]][,"current_abr",] / model$params$prob_sample[i,"hist_abr_exp"]
      
      model[["sim"]][[j]][[i]][,"bleed_mortality_ratio",] <- temp_bleed_ratio * (model$params$prob_sample[i,"bleed_mortality_ratio"] - 1) + 1
      
      #===============================================================================
      # Total Death Probability
      #===============================================================================
      
      # The SMR pertains to the death rate, so we need to multiply the death rate (not probability) by age per cycle
      
      temp_death_rate_age <- -log(1-model[["sim"]][[j]][[i]][,"death_prob_age",])
      
      # The remaining SMR is then applied to the age-death-rate, to see by how much this increases.
      
      temp_death_rate_total <- temp_death_rate_age * model[["sim"]][[j]][[i]][,"bleed_mortality_ratio",]
      
      # This adjusted death rate is converted back into a probability
      
      model[["sim"]][[j]][[i]][,"death_prob_total",] <- (1 - exp(-temp_death_rate_total))
        
      #===============================================================================
      # Probability of being alive at the start of a cycle
      #===============================================================================
      
      # The probability of being alive at the start of a given cycle is the product of all previous (1 - probabilities) of death
      # As this accumulates from cycle to cycle, we loop over the cycles (starting with cycle 2) chronologically, rather than over the patients
      
      # Under "expected_value": the alive variable takes the probability of a patient being alive at the start of a cycle
      
      if (mode == "expected_value") {
        
        for (k in 2:model$struct$ncycle) {
          model[["sim"]][[j]][[i]][k,"alive",] <- model[["sim"]][[j]][[i]][k-1,"alive",] * (1 - model[["sim"]][[j]][[i]][k-1,"death_prob_total",])
        }
        
      # Under "random": alive is set to be 0 or 1 as a random chance given the probability of living from last cycle to the current one
        
      } else if (mode == "random") {
        
        for (k in 2:model$struct$ncycle) {
          model[["sim"]][[j]][[i]][k,"alive",] <- model[["sim"]][[j]][[i]][k-1,"alive",] * rbinom(n=length(model[["sim"]][[j]][[i]][k,"alive",]), size=1, prob=(1 - model[["sim"]][[j]][[i]][k-1,"death_prob_total",]))
        }
        
      }
      
      #===============================================================================
      # Probability of being alive for the cycle as a whole
      #===============================================================================
      
      # If we base our life year calculations on the probability that they are alive at the start of it, 
      # we will overestimate survival, since patients may die at any point during the cycle
      # We correct for this by taking the average between the probability between the current and the next cycle
      
      for (k in 1:(model$struct$ncycle-1)) {
        model[["sim"]][[j]][[i]][k,"alive_corrected",] <- (model[["sim"]][[j]][[i]][k,"alive",] + model[["sim"]][[j]][[i]][k+1,"alive",]) / 2
      }
      
      #===============================================================================
      # Life Years (LYs)
      #===============================================================================
      
      # Based on probability of being alive and cycle length, we calculate undiscounted LYs accumulated per cycle
      
      model[["sim"]][[j]][[i]][,"LYs_undisc",] <- model[["sim"]][[j]][[i]][,"alive_corrected",] * (model$struct$cyc2day/settings$time$yr2day)
      
    }
  }
  
  return(model)
  
}
