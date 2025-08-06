#' Function to save selected results of models in data frame
#' Includes: 
#' 1. Outcomes
#' 2. Analysis

#' If mode is "simulation": 
#' The probabilistic sensitivity analysis can require a lot of memory if many 
#' simulations are generated with a large number of patients.
#' Therefore, it may be necessary to run some simulations, save the outputs,
#' and then delete the model from the environment before running more simulations.
#' This function lets you do this.

#' If mode is "patient": 
#' The function can be used to check convergence of mean results for patients
#' in the model. This is useful to check if the sample size is large enough
#' per simulation to get a stable result.
#' This only checks the first simulation, even if multiple ones are included in 
#' the model structure. Therefore it makes most sense to use with a deterministic
#' model.

fun_results <- function(model, mode){ 
  
  #===============================================================================
  # 1. Outcomes
  #===============================================================================
  
  # We create a vector of the outcomes we want to track
  
  outcomes <- c("cost","LYs","QALYs",
                "cost_ETRANACOGENE", "cost_PROPHYLAXIS",
                "cost_factor_treat_bleed", "cost_hospitalization", "cost_surgery")
  
  # We create an empty vector that can hold the names of all the outcomes for each of our treatment arms.
  
  intervention_outcomes <-  vector(mode="character", length = length(outcomes) * length(model[["sim"]])  )
  
  # We append the names of the treatment arms to the outcomes, so that we can save the
  # outcomes of all treatment arms in a single data frame and still tell them apart.
  
  for (i in 1:length(model[["sim"]])) {
    for (j in 1:length(outcomes)) {
      intervention_outcomes[(i-1) * length(outcomes) + j] <- paste0(names(model[["sim"]][i]),"_",outcomes[j],"_undisc")
      intervention_outcomes[(i-1) * length(outcomes) + j + length(outcomes)*length(model[["sim"]])] <- paste0(names(model[["sim"]][i]),"_",outcomes[j],"_disc")
    }
  }
  
  # We extract the outcomes we are interested in into our table
  # If the mode is "simulation": We take the mean of the individuals per simulation for each outcome
  # If the mode is "patient": We take the final value of the individuals for each outcome
  
  if (mode == "simulation") {
    
    # We construct a data frame to compare the outcomes of our probabilistic models.
    # In addition to the outcomes of our model, we will later use this to calculate ICER and Net Monetary Benefit (NMB).
    # We will also add the parameters of our probabilistic simulations later.
    # Rows: Simulations (1 per simulation)
    # Columns: Outcomes/Parameters
    
    model_outcomes <- setNames(data.frame(matrix(data = 0,
                                                 nrow = model$struct$nsim,
                                                 ncol = length(intervention_outcomes))), 
                               intervention_outcomes)
    
    # i = Number of simulations
    # j = Number of treatment arms
    # k = Number of outcomes per treatment arm
    
    # We add the undiscounted outcomes to the data frame
    
    for (i in 1:model$struct$nsim) {
      for (j in 1:length(model$sim)) {
        for (k in 1:length(outcomes))  {
          model_outcomes[i, (j-1) * length(outcomes) + k] <- sum(model[["sim"]][[j]][[i]][,paste0(outcomes[k],"_undisc"),]) / model$struct$npatients
          
        } 
      }
    }
    
    # We calculate the discounted outcomes
    
    for (i in 1:model$struct$nsim) {
      for (j in 1:length(model$sim)) {
        for (k in 1:length(outcomes))  {
          model_outcomes[i, (j-1) * length(outcomes) + k + length(outcomes)*length(model[["sim"]])] <- sum(model[["sim"]][[j]][[i]][,paste0(outcomes[k],"_undisc"),]  * model[["disc_factors"]][["disc_factors"]]) / model$struct$npatients
          
        } 
      }
    }
    
    
  } else if (mode == "patient") {
    
    # We construct a data frame to compare the outcomes of our patients.
    # In addition to the outcomes of our model, we will later use this to calculate ICER and Net Monetary Benefit (NMB).
    # Rows: Patients (1 per patient)
    # Columns: Outcomes
    
    model_outcomes <- setNames(data.frame(matrix(data = 0,
                                                 nrow = model$struct$npatients*model$struct$nsim,
                                                 ncol = length(intervention_outcomes))), 
                               intervention_outcomes)
    
    
    # i = Number of simulations
    # j = Number of treatment arms
    # k = Number of outcomes per treatment arm
    
    
    # We add the undiscounted outcomes to the data frame
    
    for (i in 1:model$struct$nsim) {
      for (j in 1:length(model$sim)) {
        for (k in 1:length(outcomes))  {
          model_outcomes[(1 + (i-1) * model$struct$npatients) : (i * model$struct$npatients), (j-1) * length(outcomes) + k] <- colSums(model[["sim"]][[j]][[i]][,paste0(outcomes[k], "_undisc"),])
            
        }
      }
    }
    
    # We calculate the discounted outcomes
    
    for (i in 1:model$struct$nsim) {
        for (j in 1:length(model$sim)) {
          for (k in 1:length(outcomes))  {
            model_outcomes[(1 + (i-1) * model$struct$npatients) : (i * model$struct$npatients), (j-1) * length(outcomes) + k + length(outcomes)*length(model[["sim"]])] <- colSums(model[["sim"]][[j]][[i]][,paste0(outcomes[k], "_undisc"),] * model[["disc_factors"]][["disc_factors"]])
            
        }
      }
    }
  }
  
  #===============================================================================
  # 2. Analysis
  #===============================================================================
  
  #' We are interested in the
  #' Net Monetary Benefit (NMB)
  #' Incremental Net Monetary Benefit (INMB)
  #' Incremental QALYs
  #' Incremental Costs
  #' Incremental Cost-Effectiveness Ratio (ICER).
  
  ## NMB
  
  # We again create a vector of names and an empty data frame to hold our NMB
  
  intervention_NMB <-  vector(mode="character", length = length(model[["sim"]]))
  
  for (i in 1:length(model[["sim"]])) {
    intervention_NMB[i] <- paste0(names(model[["sim"]][i]),"_NMB")
  }
  
  if (mode == "simulation") {
    
    model_NMB <- setNames(data.frame(matrix(data = 0,
                                            nrow = model$struct$nsim,
                                            ncol = length(intervention_NMB))),
                          intervention_NMB)
    
  } else if (mode == "patient") {
    
    model_NMB <- setNames(data.frame(matrix(data = 0,
                                            nrow = model$struct$npatients*model$struct$nsim,
                                            ncol = length(intervention_NMB))),
                          intervention_NMB)
    
  }
  
  # We calculate the Net Monetary Benefit (NMB) for each simulation in each treatment arm.
  # The formula is the following: NMB = Discounted QALYs * WTP-Threshold - Discounted Costs
  
  for (i in 1:length(model[["sim"]])) {
    model_NMB[i] <- (model_outcomes[,paste0(names(model[["sim"]][i]),"_QALYs_disc")] * model$struct$threshold) - model_outcomes[,paste0(names(model[["sim"]][i]),"_cost_disc")]
  }
  
  #' For comparative measures (INMB, incremental QALYS, incremental Costs, and ICER),
  #' later treatment are compared to earlier treatments as listed in the model structure.
  #' So if we have 4 treatments which are listed in the order X, W, Y, Z and we calculate the ICER:
  #' ICER 1: X vs W
  #' ICER 2: X vs Y
  #' ICER 3: X vs Z
  #' ICER 4: W vs Y
  #' ICER 5: W vs Z
  #' ICER 6: Y vs Z
  #' The calculation would then be: ICER 1 = (Costs X - Costs W) / (QALYs X - QALYs W)
  #' They are named in the format "X_W_ICER"
  
  ## INMB
  
  # We need a number of names equal to the number of comparisons between all the treatment arms
  # The number of comparisons (S) for the number of treatments (n) is equal to:
  # S = (n-1)*(n)/2 
  
  intervention_INMB <-  vector(mode="character", length = (length(model[["sim"]])-1)*length(model[["sim"]])/2)
  
  # We create names for the pairwise comparisons
  
  k <- 1
  
  for (i in 1:(length(model[["sim"]])-1)) {
    for (j in (i+1):length(model[["sim"]])) {
      intervention_INMB[k] <- paste0(names(model[["sim"]][i]),"_", names(model[["sim"]][j]),"_INMB")
      k <- k+1
    }
  }
  
  # Empty data frame with a column for each comparison
  
  if (mode == "simulation") {
    
    model_INMB <- setNames(data.frame(matrix(data = 0,
                                             nrow = model$struct$nsim,
                                             ncol = length(intervention_INMB))),
                           intervention_INMB)
    
    
  } else if (mode == "patient") {
    
    model_INMB <- setNames(data.frame(matrix(data = 0,
                                             nrow = model$struct$npatients*model$struct$nsim,
                                             ncol = length(intervention_INMB))),
                           intervention_INMB)
    
    
  }
  
  # We calculate the INMB
  
  k <- 1
  
  for (i in 1:(length(model[["sim"]])-1)) {
    for (j in (i+1):length(model[["sim"]])) {
      model_INMB[k] <- model_NMB[i] - model_NMB[j]
      k <- k+1
    }
  }
  
  ## Incremental QALYs, Costs, and ICER
  
  # We need a number of names equal to the number of comparisons between all the treatment arms
  # The number of comparisons (S) for the number of treatments (n) is equal to:
  # S = (n-1)*(n)/2 
  
  intervention_inc_QALY   <- vector(mode="character", length = (length(model[["sim"]])-1)*length(model[["sim"]])/2)
  intervention_inc_Costs  <- vector(mode="character", length = (length(model[["sim"]])-1)*length(model[["sim"]])/2)
  intervention_ICER       <- vector(mode="character", length = (length(model[["sim"]])-1)*length(model[["sim"]])/2)
  
  # We create names for the pairwise comparisons
  
  k <- 1
  
  for (i in 1:(length(model[["sim"]])-1)) {
    for (j in (i+1):length(model[["sim"]])) {
      
      intervention_inc_QALY[k]  <- paste0(names(model[["sim"]][i]),"_", names(model[["sim"]][j]),"_inc_QALY")
      intervention_inc_Costs[k] <- paste0(names(model[["sim"]][i]),"_", names(model[["sim"]][j]),"_inc_Costs")
      intervention_ICER[k]      <- paste0(names(model[["sim"]][i]),"_", names(model[["sim"]][j]),"_ICER")
      
      k <- k+1
    }
  }
  
  # Empty data frame with a column for each comparison
  
  if (mode == "simulation") {
    
    model_inc_QALY  <- setNames(data.frame(matrix(data = 0,
                                                  nrow = model$struct$nsim,
                                                  ncol = length(intervention_inc_QALY))),
                                intervention_inc_QALY)
    
    model_inc_Costs <- setNames(data.frame(matrix(data = 0,
                                                  nrow = model$struct$nsim,
                                                  ncol = length(intervention_inc_Costs))),
                                intervention_inc_Costs)
    
    model_ICER      <- setNames(data.frame(matrix(data = 0,
                                                  nrow = model$struct$nsim,
                                                  ncol = length(intervention_ICER))),
                                intervention_ICER)
    
  } else if (mode == "patient") {
    
    model_inc_QALY  <- setNames(data.frame(matrix(data = 0,
                                                  nrow = model$struct$npatients*model$struct$nsim,
                                                  ncol = length(intervention_inc_QALY))),
                                intervention_inc_QALY)
    
    model_inc_Costs <- setNames(data.frame(matrix(data = 0,
                                                  nrow = model$struct$npatients*model$struct$nsim,
                                                  ncol = length(intervention_inc_Costs))),
                                intervention_inc_Costs)
    
    model_ICER      <- setNames(data.frame(matrix(data = 0,
                                                  nrow = model$struct$npatients*model$struct$nsim,
                                                  ncol = length(intervention_ICER))),
                                intervention_ICER)
    
  }
  
  # We calculate the incremental QALYs, Costs, and ICER
  
  k <- 1
  
  for (i in 1:(length(model[["sim"]])-1)) {
    for (j in (i+1):length(model[["sim"]])) {
      
      model_inc_QALY[k]  <- model_outcomes[,paste0(names(model[["sim"]][i]),"_QALYs_disc")] - model_outcomes[,paste0(names(model[["sim"]][j]),"_QALYs_disc")]
      model_inc_Costs[k] <- model_outcomes[,paste0(names(model[["sim"]][i]),"_cost_disc")] - model_outcomes[,paste0(names(model[["sim"]][j]),"_cost_disc")]
      model_ICER[k]      <- model_inc_Costs[k] / model_inc_QALY[k]
      
      k <- k+1
    }
  }
  
  # We combine all of our analyses
  
  model_analysis <- cbind(model_NMB, model_INMB, model_inc_QALY, model_inc_Costs, model_ICER)
  
  
  #===============================================================================
  # 3. Patient characteristics
  #===============================================================================
  
  # If we are looking at individual patients, extract most relevant characteristics
  # of these patients.
  # These are: 
  # 1. Age start of model
  # 2. Sex
  # 3. Baseline ABR (individual)
  # As these characteristics are independent of treatment assignment, only
  # need to do this once per treatment.
  
  if (mode == "patient") {

    model_patient_char <- setNames(data.frame(matrix(data = 0,
                                                 nrow = model$struct$npatients*model$struct$nsim,
                                                 ncol = 3)),
                               c("sex", "age","baseline_abr_ind"))


      for (i in 1:model$struct$nsim) {

        model_patient_char[(1 + (i-1) * model$struct$npatients) : (i * model$struct$npatients), "sex"] <- model[["sim"]][[1]][[i]][1,"sex",]
        model_patient_char[(1 + (i-1) * model$struct$npatients) : (i * model$struct$npatients), "age"] <- model[["sim"]][[1]][[i]][1,"age",]
        model_patient_char[(1 + (i-1) * model$struct$npatients) : (i * model$struct$npatients), "baseline_abr_ind"] <- model[["sim"]][[1]][[i]][1,"baseline_abr",]

      }

  }
  
  
  #===============================================================================
  # Combine and Return
  #===============================================================================
  
  if (mode == "simulation") {
    
    # We combine outcomes, analysis, and probabilistic parameters into a single data frame
    model_results <- cbind(model_outcomes, model_analysis, model$params$prob_sample)
    
  } else if (mode == "patient") {
    
    # We assign probabilistic parameters for their simulation to each patient.
    
    param_names <- c(names(data.frame(model$params$prob_sample)))
    
    simulation_parameters <- setNames(data.frame(matrix(data = 0,
                                                        nrow = model$struct$npatients*model$struct$nsim,
                                                        ncol = length(param_names))), 
                                      param_names)
    
    for (i in 1:model$struct$nsim) {
      for (j in 1:model$struct$npatients)
        
        simulation_parameters[(i-1)*model$struct$npatients+j,] <- model$params$prob_sample[i,]
      
    }
    
    # We combine outcomes, patient characteristics, analysis, and probabilistic parameters into a single data frame
    model_results <- cbind(model_outcomes, model_patient_char, model_analysis, simulation_parameters)

    
  }
  
  # We return our combined results from the function
  return(model_results)
  
}