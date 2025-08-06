#' Function to determine treatment and number of bleeds per cycle in the ONDEMAND treatment arm
#' Treatment 31: On-Demand treatment

fun_ONDEMAND_abr <- function(model){ 
  
  for (i in 1:model$struct$nsim) {
    
    # Apply treatment
    model[["sim"]][["ONDEMAND"]][[i]][,"treatment", ]  <- 31
    
    # Apply no bleed reduction in all periods for all patients in treatment arm
    model[["sim"]][["ONDEMAND"]][[i]][,"abr_reduction",] <- 0
    
    # Receive ABR by multiplying baseline ABR with bleed reduction in each period
    model[["sim"]][["ONDEMAND"]][[i]][,"current_abr", ] <- model[["sim"]][["ONDEMAND"]][[i]][,"baseline_abr", ] * (1 - model[["sim"]][["ONDEMAND"]][[i]][,"abr_reduction", ])
    
  }
  
  return(model)
  
}
