#' Function to set scenarios for the model
#' List of Scenarios:
#' 1. Treatment Switching
#' 2. Vial Sharing

fun_scenarios <- function(model){
  
  ################################################################################
  #                                                                              #
  #  SCENARIOS                                                                   #
  #                                                                              #
  ################################################################################
  
  # Set binary (0/1) switches to determine whether scenario applies or not
  
  #-------------------------------------------------------------------------------
  # Treatment Switching
  #-------------------------------------------------------------------------------
  
  # If sc_treatment_switching = 1, then patients receiving Etranacogene Dezaparvovec
  # switch to prophylaxis if current_abr exceed model$scenarios$treatment_switch_abr
  
  model$scenarios$sc_treatment_switching  <- 1
  
  model$scenarios$treatment_switch_abr    <- 4
  
  #-------------------------------------------------------------------------------
  # Vial Sharing
  #-------------------------------------------------------------------------------
  
  # If sc_vial_sharing = 1, then patients share vials, meaning none of the coagulation factor is wasted
  # This reduces the resource use for PROPHYLAXIS and for infusions to treat bleeds, making treatment more efficient and thus cheaper
  # This assumption is realistic if patients are treated in centers, and probably less realistic for home treatment
  
  model$scenarios$sc_vial_sharing  <- 1
  
  # If patients do not share vials, then vial size is important. 
  # Larger vial size leads to more waste, as we round up to next largest vial size.
  # See EMA 2017: https://www.ema.europa.eu/en/documents/product-information/refixia-epar-product-information_en.pdf
  # Vial size 1: 500 IU (smallest)
  # Vial size 2: 1000 IU
  # Vial size 3: 2000 IU (largest)
  # By default, we assume the smallest vial size of 500 IU
  
  model$scenarios$vial_size <- 500
  
  #-------------------------------------------------------------------------------
  # Longer bleed disutility
  #-------------------------------------------------------------------------------
  
  # In the base case, we assume a bleed disutility duration of 1 day per bleed Neufeld et al. (2012)
  # This may not be entirely realistic, so we test a longer disutility duration in a scenario
  # By default, we assume the utility lasts one day
  
  model$scenarios$sc_bleed_disutility_one_day  <- 1
  
  # In the scenario, we test the effects of it lasting 7 days.
  
  model$scenarios$bleed_disutility_duration <- 7
  
  ################################################################################
  #                                                                              #
  #  END                                                                         #
  #                                                                              #
  ################################################################################
  
  return(model)
  
}