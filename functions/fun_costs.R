#' Function to determine costs based on resource use and prices

fun_costs <- function(model){ 
  
  for (i in 1:model$struct$nsim) {
    for (j in 1:length(model$sim)) {
      
      #===============================================================================
      # Treatment
      #===============================================================================
      
      # Multiply treatment with viral gene Therapy by price of treatment
      model[["sim"]][[j]][[i]][,"cost_ETRANACOGENE_undisc",] <- model[["sim"]][[j]][[i]][,"res_ETRANACOGENE",] * model[["params"]][["prob_sample"]][i,"price_ETRANACOGENE"]
      
      # Multiply treatment with Prophylaxis (based on quantity of factor) by price of coagulation factor
      model[["sim"]][[j]][[i]][,"cost_PROPHYLAXIS_undisc",] <- model[["sim"]][[j]][[i]][,"res_IU_PROPHYLAXIS",] * model[["params"]][["prob_sample"]][i,"price_coagulation_factor_IU"]
      
      #===============================================================================
      # Coagulation Factor to treat bleeds
      #===============================================================================
      
      # Multiply treatment with coagulation factor to treat bleeds (based on quantity of factor) by price of coagulation factor
      model[["sim"]][[j]][[i]][,"cost_factor_treat_bleed_undisc",] <- model[["sim"]][[j]][[i]][,"res_IU_treat_bleed",] * model[["params"]][["prob_sample"]][i,"price_coagulation_factor_IU"]
      
      #===============================================================================
      # Hospitalization
      #===============================================================================
      
      # Multiply number of hospitalizations by price per hospitalization, summing up ward and ICU days
      model[["sim"]][[j]][[i]][,"cost_hospitalization_undisc",] <- model[["sim"]][[j]][[i]][,"res_hosp_ward_days",] * model[["params"]][["prob_sample"]][i,"price_hosp_ward_days"] + model[["sim"]][[j]][[i]][,"res_hosp_ICU_days",] * model[["params"]][["prob_sample"]][i,"price_hosp_ICU_days"]
      
      #===============================================================================
      # Surgery
      #===============================================================================
      
      # Multiply number of surgeries by price per surgery
      # We also assume that a surgery requires an infusion of FIX and add these costs as well
      model[["sim"]][[j]][[i]][,"cost_surgery_undisc",] <- model[["sim"]][[j]][[i]][,"res_surgery",] * model[["params"]][["prob_sample"]][i,"price_surgery"] + model[["sim"]][[j]][[i]][,"res_surgery",] * model[["sim"]][[j]][[i]][,"weight",] * model$params$prob_sample[i,"res_IU_kg_joint_bleed"] * model[["params"]][["prob_sample"]][i,"price_coagulation_factor_IU"]
      
      #===============================================================================
      # Total Costs
      #===============================================================================
      
      # We sum up all the other undiscounted costs
      model[["sim"]][[j]][[i]][,"cost_undisc",] <- model[["sim"]][[j]][[i]][,"cost_ETRANACOGENE_undisc",] + model[["sim"]][[j]][[i]][,"cost_PROPHYLAXIS_undisc",] + model[["sim"]][[j]][[i]][,"cost_factor_treat_bleed_undisc",] + model[["sim"]][[j]][[i]][,"cost_hospitalization_undisc",] + model[["sim"]][[j]][[i]][,"cost_surgery_undisc",] 
      
    }
  }
  
  return(model)
  
}