# Function to return tornado plot of univariate sensitivity analysis
# Requires 3 inputs:
# 1. Sensitivity analysis list
# 2. outcome of choice for tornado plot
# 3. Baseline value of that outcome

fun_tornado <- function(sens_list, outcome, baseline){ 
  
  # Determine parameters to show in tornado diagram
  
  parameters <- names(sens_analysis)
  
  # Determine baseline value of outcome
  
  baseline <- baseline_model_results[,outcome]
  
  # Create empty data frame to hold data with required columns
  
  tornado_data <- data.frame(matrix(nrow = length(parameters) * 2, ncol = 11))
  colnames(tornado_data) <- c("parameter", "type", "outcome_value", "parameter_value", "difference","ymin","ymax","xmin","xmax", "label_pos", "rank")
  
  # Fill tornado with data
  
  for (i in 1:length(parameters)) {
    
    # Assign two lines to each parameter in the sensitivity analysis, one for high and one for low
    
    tornado_data[1+((i-1)*2),"parameter"] <- parameters[i]
    tornado_data[2+((i-1)*2),"parameter"] <- parameters[i]
    
    # Designate line for upper bound and find highest value of that parameter in the sensitivity analysis
    
    tornado_data[1+((i-1)*2),"type"] <- "upper_bound"
    tornado_data[1+((i-1)*2),"parameter_value"] <- max(sens_analysis[[paste0(parameters[i])]][["range"]])
    
    # Determine position of that parameter value and assign parameter value to data set
    
    pos_upper <- which(tornado_data[1+((i-1)*2),"parameter_value"] == sens_analysis[[paste0(parameters[i])]][["range"]])
    tornado_data[1+((i-1)*2),"outcome_value"] <- sens_analysis[[paste0(parameters[i])]][["model_results"]][paste0(outcome)][pos_upper,]
    
    # Designate line for lower bound and find lowest value of that outcome in the sensitivity analysis
    
    tornado_data[2+((i-1)*2),"type"] <- "lower_bound"
    tornado_data[2+((i-1)*2),"parameter_value"] <- min(sens_analysis[[paste0(parameters[i])]][["range"]])
    
    # Determine position of that parameter value and assign parameter value to data set
    
    pos_lower <- which(tornado_data[2+((i-1)*2),"parameter_value"] ==sens_analysis[[paste0(parameters[i])]][["range"]])
    tornado_data[2+((i-1)*2),"outcome_value"] <-  sens_analysis[[paste0(parameters[i])]][["model_results"]][paste0(outcome)][pos_lower,]
    
    # Determine difference of outcome between upper and lower
    
    tornado_data[1+((i-1)*2),"difference"]  <- max(tornado_data[1+((i-1)*2),"outcome_value"],tornado_data[2+((i-1)*2),"outcome_value"]) - min(tornado_data[1+((i-1)*2),"outcome_value"], tornado_data[2+((i-1)*2),"outcome_value"])
    tornado_data[2+((i-1)*2),"difference"]  <- tornado_data[1+((i-1)*2),"difference"]
    
    # Set ymin and ymax for upper bound
    
    tornado_data[1+((i-1)*2),"ymin"] <- pmin(tornado_data[1+((i-1)*2),"outcome_value"], baseline)
    tornado_data[1+((i-1)*2),"ymax"] <- pmax(tornado_data[1+((i-1)*2),"outcome_value"], baseline)
    
    # Set ymin and ymax for lower bound
    
    tornado_data[2+((i-1)*2),"ymin"] <- pmin(tornado_data[2+((i-1)*2),"outcome_value"], baseline)
    tornado_data[2+((i-1)*2),"ymax"] <- pmax(tornado_data[2+((i-1)*2),"outcome_value"], baseline)
    
    # Determine position of upper bound label
    
    if (tornado_data[1+((i-1)*2),"outcome_value"] > tornado_data[2+((i-1)*2),"outcome_value"]){
  
      tornado_data[1+((i-1)*2),"rank"] <- "high"
    } else {
      tornado_data[1+((i-1)*2),"label_pos"] <- tornado_data[1+((i-1)*2),"ymax"]
      tornado_data[1+((i-1)*2),"rank"] <- "low"
    }
    
    # Determine position of lower bound label
    
    if (tornado_data[2+((i-1)*2),"outcome_value"] < tornado_data[1+((i-1)*2),"outcome_value"]){
      tornado_data[2+((i-1)*2),"label_pos"] <- tornado_data[2+((i-1)*2),"ymax"]
      tornado_data[2+((i-1)*2),"rank"] <- "low"
    } else {
      tornado_data[2+((i-1)*2),"label_pos"] <- tornado_data[2+((i-1)*2),"ymin"]
      tornado_data[2+((i-1)*2),"rank"] <- "high"
    }
    
    tornado_data[1+((i-1)*2),"label_pos"] <- tornado_data[1+((i-1)*2),"ymin"]
  }
  
  # Sort dataset by size of difference and get ordered vector of parameter names
  
  tornado_data     <- tornado_data[order(tornado_data$type, tornado_data$difference, tornado_data$parameter),]
  order_parameters <- tornado_data$parameter[order(tornado_data$type, tornado_data$difference, tornado_data$parameter)]
  
  # Remove duplicate parameter names
  
  for (i in 1:length(parameters)) {
    order_parameters <- order_parameters[-2]
  }
  
  # Width of columns in plot (value between 0 and 1)
  width <- 0.95
  
  # Set xmin and xmax depending on order in dataset
  
  for (i in 1:length(parameters)) {
    
    # Set xmin and xmax for lower bound
    
    tornado_data[i,"xmin"] <- i - width/2
    tornado_data[i,"xmax"] <- i + width/2
    
    # Set xmin and xmax for upper bound
    
    tornado_data[i+length(order_parameters),"xmin"] <- i - width/2
    tornado_data[i+length(order_parameters),"xmax"] <- i + width/2
    
  }
  
  return(tornado_data)
  
}

