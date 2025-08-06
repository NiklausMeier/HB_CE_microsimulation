#' Function to generate decision rules data structure
#' and determine optimal treatment assignment based on counterfactuals
#' Input: 
#' 1. Data from microsimulation
#' 2. Number of treated patients
#' 3. Vector of treatments
#' 4. Patient characteristics to include as covariates
#' 5. Outcome that should be optimized for

fun_opt_treat <- function(data, ntreat, treatments, covariates, optimize){ 

  #===============================================================================
  # Prepare dec_rules list
  #===============================================================================
  
  # Create list to hold relevant data
  dec_rules <- list()
  
  # Are we trying to optimize costs, QALYs, or NMB?
  dec_rules$optimize <- optimize
  
  # Which patient characteristics are we using to optimize treatment?
  dec_rules$covariates <- covariates
  
  # Copy over patient results
  dec_rules$PRC <- data
  
  # Number of treatments
  dec_rules$ntreat <- ntreat
  
  # n = number of observed patients. We multiply the original number of patients by
  # the number of treatments since we treat each patient as separate when 
  # "randomized" into the different treatments.
  dec_rules$n <- nrow(PRC) * dec_rules$ntreat
  
  # Define our covariates X in a matrix
  dec_rules$X <- matrix(0, dec_rules$n, 3)
  colnames(dec_rules$X) <- dec_rules$covariates
  
  # Outcomes as a vector
  dec_rules$Y <- vector(length = dec_rules$n)
  
  # We create a factor W with levels for each of the treatments.
  dec_rules$W <- factor(levels=treatments)
  
  # We generate bins of observations for each of our treatments to help extract
  # covariates and outcomes from our data.
  dec_rules$bins <- vector("list", length(model[["sim"]]))
  
  # we loop over the number of treatments to extract treatments (W),
  # covariates (X), and outcomes (Y) from the patient results
  for (i in 1:dec_rules$ntreat) {
    
    # Define bin i
    dec_rules$bins[[i]] <- seq(from = ((i-1)*nrow(dec_rules$PRC) + 1), to = (i * nrow(dec_rules$PRC)), by = 1)
    
    # Treatments for patients in bin i
    dec_rules$W[dec_rules$bins[[i]]]    <- levels(dec_rules$W)[i]
    
    # Covariates for patients in bin i
    dec_rules$X[dec_rules$bins[[i]],"baseline_abr_ind"] <- round(dec_rules$PRC[,"baseline_abr_ind"],0)
    dec_rules$X[dec_rules$bins[[i]],"age"]              <- round(dec_rules$PRC[,"age"],0)
    dec_rules$X[dec_rules$bins[[i]],"sex"]              <- dec_rules$PRC[,"sex"]
    
    # Outcomes Y for patients in bin i
    if (optimize == "NMB") {
      
      dec_rules$Y[dec_rules$bins[[i]]]   <- round(dec_rules$PRC[,paste0(levels(dec_rules$W)[i],"_NMB")],0)
      
    } else if (optimize == "QALY") {
      
      dec_rules$Y[dec_rules$bins[[i]]]   <- round(dec_rules$PRC[,paste0(levels(dec_rules$W)[i],"_QALYs_disc")],2)
      
    } else if (optimize == "COST") {
      
      dec_rules$Y[dec_rules$bins[[i]]]   <- -round(dec_rules$PRC[,paste0(levels(dec_rules$W)[i],"_cost_disc")],0)
      
    }
  }
  
  #-------------------------------------------------------------------------------
  # QALY
  #-------------------------------------------------------------------------------
  
  # Create empty data frame to hold QALY data
  dec_rules$patient_QALY <- setNames(data.frame(matrix(data = 0,
                                                       nrow = nrow(dec_rules$PRC),
                                                       ncol = dec_rules$ntreat+1)), 
                                     c(levels(dec_rules$W), "OPTIMAL"))
  
  
  # Extract QALYs per patient depending on assigned treatment
  for (i in 1:dec_rules$ntreat) {
    
    dec_rules$patient_QALY[,paste0(levels(dec_rules$W)[i])] <- dec_rules$PRC[,paste0(levels(dec_rules$W)[i],"_QALYs_disc")]
    
  }
  
  #-------------------------------------------------------------------------------
  # Costs
  #-------------------------------------------------------------------------------
  
  # Create empty data frame to hold costdata
  dec_rules$patient_COST <- setNames(data.frame(matrix(data = 0,
                                                       nrow = nrow(dec_rules$PRC),
                                                       ncol = dec_rules$ntreat+1)), 
                                     c(levels(dec_rules$W), "OPTIMAL"))
  
  
  # Extract cost per patient depending on assigned treatment
  for (i in 1:dec_rules$ntreat) {
    
    dec_rules$patient_COST[,paste0(levels(dec_rules$W)[i])] <- dec_rules$PRC[,paste0(levels(dec_rules$W)[i],"_cost_disc")]
    
  }
  
  #-------------------------------------------------------------------------------
  # NMB
  #-------------------------------------------------------------------------------
  
  # Create empty data frame to hold NMB data
  dec_rules$patient_NMB <- setNames(data.frame(matrix(data = 0,
                                                      nrow = nrow(dec_rules$PRC),
                                                      ncol = dec_rules$ntreat+1)), 
                                    c(levels(dec_rules$W), "OPTIMAL"))
  
  
  # Extract NMB per patient depending on assigned treatment
  for (i in 1:dec_rules$ntreat) {
    
    dec_rules$patient_NMB[,paste0(levels(dec_rules$W)[i])] <- dec_rules$PRC[,paste0(levels(dec_rules$W)[i],"_NMB")]
    
  }
  
  #===============================================================================
  # Optimal Treatment
  #===============================================================================
  
  # Assign optimal treatment based on optimization goal (QALYs, costs, or NMB)
  # From this, derive calculate QALYs, costs and NMB for each patient based on assigned treatment
  
  # Depending on what should be optimized, extract correct outcome
  
  if (optimize == "QALY") {
    
    opt_outcome <- dec_rules$patient_QALY
    
  } else if (optimize == "COST") {
    
    opt_outcome <- dec_rules$patient_COST
    
  } else if (optimize == "NMB") {
    
    opt_outcome <- dec_rules$patient_NMB
    
  }
  
  # Create empty list for pmax/pmin function to select most optimal treatment
  opt_list <- list()
  
  # Create list of vectors from data frame
  for (i in 1:dec_rules$ntreat) {
    opt_list[[i]] <- opt_outcome[,i]
  }
  
  # Pass list of arguments into pmax/pmin function via do.call
  if (optimize == "QALY") {
    
    dec_rules$patient_QALY[,"OPTIMAL"] <- do.call(pmax, opt_list)
    
  } else if (optimize == "COST") {
    
    dec_rules$patient_COST[,"OPTIMAL"] <- do.call(pmin, opt_list)
    
  } else if (optimize == "NMB") {
    
    dec_rules$patient_NMB[,"OPTIMAL"] <- do.call(pmax, opt_list)
    
  }
  
  # Create empty vector to hold treatment decisions
  dec_rules$opt_treat <- setNames(data.frame(matrix(data = 0,
                                                    nrow = nrow(dec_rules$PRC),
                                                    ncol = dec_rules$ntreat+1)), 
                                  c("OPTIMAL", levels(dec_rules$W)))
  
  # Check which treatment this corresponds to and assign name of optimal treatment
  for (i in 1:dec_rules$ntreat) {
    
    
    if (optimize == "QALY") {
      
      treat <- which(dec_rules$patient_QALY[,i] == dec_rules$patient_QALY[,"OPTIMAL"])
      
    } else if (optimize == "COST") {
      
      treat <- which(dec_rules$patient_COST[,i] == dec_rules$patient_COST[,"OPTIMAL"])
      
    } else if (optimize == "NMB") {
      
      treat <- which(dec_rules$patient_NMB[,i] == dec_rules$patient_NMB[,"OPTIMAL"])
      
    }
    
    dec_rules$opt_treat[treat,"OPTIMAL"] <- levels(dec_rules$W)[i]
    dec_rules$opt_treat[treat,i+1] <- 1
    
  }
  
  psum2 <- function(...,na.rm=FALSE) { 
    dat <- do.call(cbind,list(...))
    res <- rowSums(dat, na.rm=na.rm) 
    idx_na <- !rowSums(!is.na(dat))
    res[idx_na] <- NA
    res 
  }
  
  # For QALYs, costs, and NMB, assign treatment based on optimal treatment assignment
  dec_rules[["patient_QALY"]][,"OPTIMAL"] <- psum2(dec_rules$opt_treat[,-1] * dec_rules[["patient_QALY"]][,1:dec_rules$ntreat])
  dec_rules[["patient_COST"]][,"OPTIMAL"] <- psum2(dec_rules$opt_treat[,-1] * dec_rules[["patient_COST"]][,1:dec_rules$ntreat])
  dec_rules[["patient_NMB"]][,"OPTIMAL"]  <- psum2(dec_rules$opt_treat[,-1] * dec_rules[["patient_NMB"]][,1:dec_rules$ntreat])
  
  #===============================================================================
  # Comparison Table
  #===============================================================================
  
  dec_rules$comp_table <- setNames(data.frame(matrix(data = 0,
                                       nrow = 5,
                                       ncol = 3)),
                                    c("QALYs", "Costs", "NMB"))
  
  rownames(dec_rules$comp_table) <- c(levels(dec_rules$W), "OPTIMAL")
  
  for (i in 1:nrow(dec_rules$comp_table)) {
    
    dec_rules$comp_table[i,"QALYs"] <- mean(dec_rules[["patient_QALY"]][[paste0(rownames(dec_rules$comp_table)[i])]])
    dec_rules$comp_table[i,"Costs"] <- mean(dec_rules[["patient_COST"]][[paste0(rownames(dec_rules$comp_table)[i])]])
    dec_rules$comp_table[i,"NMB"]   <- mean(dec_rules[["patient_NMB"]][[paste0(rownames(dec_rules$comp_table)[i])]]) 
    
  }
  
  # Determine which is the optimal uniform treatment strategy
  dec_rules$opt_uniform <- row.names(dec_rules$comp_table)[which.max(dec_rules$comp_table[1:dec_rules$ntreat,"NMB"])]
  
  #===============================================================================
  # End
  #===============================================================================
  
  # Return list with output
  return(dec_rules)
  
}
  