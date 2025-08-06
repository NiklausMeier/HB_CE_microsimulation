#' Function to determine treatment and number of bleeds per cycle in the PROPHYLAXIS treatment arm
#' Treatment 21: Prophylaxis treatment

fun_PROPHYLAXIS_abr <- function(model){ 
  
  for (i in 1:model$struct$nsim) {
    
    # Apply treatment
    model[["sim"]][["PROPHYLAXIS"]][[i]][,"treatment", ]  <- 21
    
    # Apply bleed reduction from Prophylaxis in all periods for all patients in treatment arm
    model[["sim"]][["PROPHYLAXIS"]][[i]][,"abr_reduction",] <- model$params$prob_sample[i,"PROPHYLAXIS_relative_bleed_reduction"]
    
    # Receive ABR by multiplying baseline ABR with bleed reduction in each period
    model[["sim"]][["PROPHYLAXIS"]][[i]][,"current_abr", ] <- model[["sim"]][["PROPHYLAXIS"]][[i]][,"baseline_abr", ] * (1 - model[["sim"]][["PROPHYLAXIS"]][[i]][,"abr_reduction", ])
    
  }
  
  return(model)
  
}
