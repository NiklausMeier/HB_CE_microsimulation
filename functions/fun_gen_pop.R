#' Function to generate a population for all treatment arms
#' Input: 
#' 1.   Requires a model structure that has been generated with "fun_model_setup"
#' 2.   Requires a mode to be set for population as "heterogeneous" or "homogeneous"
#' 2.A. Under "heterogeneous": Every individual has their characteristics randomly drawn
#' 2.B. Under "homogeneous": Every individual is the same, and has the expected value of the characteristics

#' The "homogeneous" option can be used to see if the model predictions differ from the average individual
#' It is essentially as though we turn our microsimulation into a cohort model

fun_gen_pop <- function(model, mode){ 
  
  # We generate random populations based on their sex, age, cumulative joint bleeds, and baseline annualized bleed rate (pre-treatment)
  # We then assign these populations to one simulation per treatment arm, so that they share populations
  # Example: If we are evaluating 10 simulations per treatment arm, we would generate 10 populations
  # Population 1 is assigned to simulation 1 of every treatment arm, population 2 is assigned to simulation 2 of every treatment arm, etc.
  
  for (i in 1:model$struct$nsim) {
    
    randomdraw <- list()
    
    #===============================================================================
    # Life
    #===============================================================================
    
    # Patients are alive at start of model
    
    for (j in 1:length(model$sim)) {
      model[["sim"]][[j]][[i]][1,"alive",]                 <- 1
    }
    
    #===============================================================================
    # Sex
    #===============================================================================
    
    if (mode == "heterogeneous") {
      
      # Draw sex of each patient
      # Female = 0
      # Male = 1
      
      randomdraw$sex                 <- sample(c(0,1),
                                               model$struct$npatients, 
                                               TRUE, 
                                               c(model$params$prob_sample[i,"share_female"],
                                                 1-model$params$prob_sample[i,"share_female"]))
      
    } else if (mode == "homogeneous") {
      
      # Almost all hemophilia patients are male (1)
      
      randomdraw$sex <- rep(1,model$struct$npatients)
      
    }
    
    # Assign sex to patients, across all cycles
    
    for (j in 1:length(model$sim)) {
      for (k in 1:model$struct$npatients) {
        model[["sim"]][[j]][[i]][,"sex",k]           <- randomdraw$sex[[k]]
      }
    }
    
    #===============================================================================
    # Age
    #===============================================================================
    
    # Create age array
    
    randomdraw$age <- array(0,dim = c(model$struct$ncycle,model$struct$npatients))
    
    
    if (mode == "heterogeneous") {
      
      # Draw age of each patient at start of simulation, in first cycle
      # Minimum age is 18, maximum age is 110
      # As we do not have an age distribution for hemophilia which is stratified by sex,
      # we assume an equal age distribution for both sexes.
      
      randomdraw$age[1,]    <- rtruncnorm(model$struct$npatients, 
                                          a = 18, 
                                          b = 110, 
                                          mean = model$params$age$mean, 
                                          sd = model$params$age$sd)
      
    } else if (mode == "homogeneous") {
      
      # Set starting age to expected value
      
      randomdraw$age[1,]    <-  model$params$age$mean       
      
    }
    
    # Increase age of patients by cycle length per cycle
    
    for (j in 1:model$struct$npatients) {
      randomdraw$age[,j]  <- randomdraw$age[1,j] + ((1:model$struct$ncycle - 1) * (model$struct$cyc2day/settings$time$yr2day))
    }
    
    # Assign age to patients
    
    for (j in 1:length(model$sim)) {
      model[["sim"]][[j]][[i]][,"age",]           <- randomdraw$age
    }
    
    #===============================================================================
    # Weight
    #===============================================================================
    
    # Create weight array
    
    randomdraw$weight <- array(0,dim = c(model$struct$ncycle,model$struct$npatients))
    
    # Determine weight based on age and sex
    
    for (j in 1:model$struct$npatients) {
      if (randomdraw$sex[[j]] == 1) {
        randomdraw$weight[,j] <- lookup(pmin(floor(randomdraw$age[,j]),110), 
                                        model$life_tables_weight$age, 
                                        model$life_tables_weight$male_weight_kg)
      } else if (randomdraw$sex[[j]] == 0) {
        randomdraw$weight[,j] <- lookup(pmin(floor(randomdraw$age[,j]),110), 
                                        model$life_tables_weight$age, 
                                        model$life_tables_weight$female_weight_kg)
      }
    }
    
    # Assign weight to patients by age and sex in each cycle
    
    for (j in 1:length(model$sim)) {
      model[["sim"]][[j]][[i]][,"weight",]          <- randomdraw$weight
    }
    
    #===============================================================================
    # Mortality based on Age
    #===============================================================================
    
    # Create death-rate-age Array
    
    randomdraw$death_rate_age <- array(0,dim = c(model$struct$ncycle,model$struct$npatients))
    
    # Determine rate of death based on age and sex
    # We do this inside the population-generation function since it needs to be appropriate to the age of each patient in the model.
    
    for (j in 1:model$struct$npatients) {
      if (randomdraw$sex[[j]] == 1) {
        randomdraw$death_rate_age[,j] <- lookup(pmin(floor(randomdraw$age[,j]), 110), 
                                                model$life_tables_weight$age, 
                                                model$life_tables_weight$male_yearly_death_rate)
      } else if (randomdraw$sex[[j]] == 0) {
        randomdraw$death_rate_age[,j] <- lookup(pmin(floor(randomdraw$age[,j]), 110), 
                                                model$life_tables_weight$age, 
                                                model$life_tables_weight$female_yearly_death_rate)
      }
    }
    
    # Create death-probability-age Array
    
    randomdraw$death_prob_age <- array(0,dim = c(model$struct$ncycle,model$struct$npatients))
    
    # Convert rate into probability, adjusted for model cycle length
    
    randomdraw$death_prob_age <- (1 - exp(-randomdraw$death_rate_age * (model$struct$cyc2day/settings$time$yr2day)))
    
    # Assign probability of death based on age to patients in each cycle
    
    for (j in 1:length(model$sim)) {
      model[["sim"]][[j]][[i]][,"death_prob_age",]          <- randomdraw$death_prob_age
    }
    
    #===============================================================================
    # Baseline annualized bleed rate (ABR)
    #===============================================================================
    
    if (mode == "heterogeneous") {
      
      # Based on our probabilistic sample mean, draw alpha and beta of gamma distribution 
      
      randomdraw$baseline_abr$alpha  <- (model$params$prob_sample[i,"baseline_abr_exp"]/model$params$baseline_abr$sd)^2
      randomdraw$baseline_abr$beta   <- (model$params$baseline_abr$sd^2)/model$params$prob_sample[i,"baseline_abr_exp"]
      
      
      # Draw random baseline annualized bleed rate
      randomdraw$baseline_abr$patients            <- rgamma(n     = model$struct$npatients, 
                                                            shape = randomdraw$baseline_abr$alpha, 
                                                            scale = randomdraw$baseline_abr$beta)
      
    } else if (mode == "homogeneous") {
      
      # Assign baseline expected ABR to all patients
      randomdraw$baseline_abr$patients <- rep(model$params$prob_sample[i,"baseline_abr_exp"], model$struct$npatients) 
      
    }
    
    # Assign baseline ABR to patients, across all cycles
    
    for (j in 1:length(model$sim)) {
      for (k in 1:model$struct$npatients) {
        model[["sim"]][[j]][[i]][,"baseline_abr",k] <- randomdraw$baseline_abr$patients[[k]]
      }
    }
    
    #===============================================================================
    # Historic annualized bleed rate (ABR)
    #===============================================================================
    
    # We get the historic ABR of all patients by multiplying their baseline ABR with our Historic/Baseline ABR Ratio
    randomdraw$hist_abr$patients                 <- randomdraw$baseline_abr$patients * model$params$prob_sample[i,"hist_baseline_abr_ratio"]
    
    #===============================================================================
    # Cumulative joint bleeds
    #===============================================================================
    
    # Calculate number of cumulative joint bleeds which a patient has at start of simulation, 
    # based on historic ABR, age, and share of joint bleeds
    
    randomdraw$bleeds_joint_cum    <- randomdraw$age[1,] * randomdraw$hist_abr$patients * model$params$prob_sample[i,"share_joint_bleeds"]
    
    # Assign number of cumulative joint bleeds to patients in first cycle
    
    for (j in 1:length(model$sim)) {
      model[["sim"]][[j]][[i]][1,"bleeds_joint_cum",]  <- randomdraw$bleeds_joint_cum
    }
    
    
    
    
  }
  
  return(model)
  
}
