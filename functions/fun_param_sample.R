#' Function to sample parameters for each simulation
#' Even if the model is deterministic, this function still needs to be run, so that deterministic parameters are drawn
#' Inputs: Model (with all the required probabilistic parameters)
#' If the model is deterministic, the mean values are applied all simulations
#' If the model is probabilistic, random values are drawn based on the parameters and distribution for each simulation.

fun_param_sample <- function(model){ 
  
  #===============================================================================
  # Prepare vector of probabilistic parameter names
  #===============================================================================
  
  # We make a vector of all probabilistic parameters we sample once per simulation
  
  model$params$prob_names <- c("util_baseline", "util_male", "util_age", "util_agesq",
                               "util_PS_13_to_21", "util_PS_22_plus", 
                               "util_joint_bleed", "util_other_bleed", 
                               "util_infusion", "util_surgery",
                               
                               "share_female",
                               "hist_abr_exp", "baseline_abr_exp", "hist_baseline_abr_ratio",
                               
                               "PS_bleeds",
                               "share_joint_bleeds", "share_other_bleeds",
                               "bleed_mortality_ratio",
                               
                               "price_ETRANACOGENE", 
                               "price_coagulation_factor_IU", "price_hosp_ward_days",
                               "price_hosp_ICU_days", "price_surgery",
                               
                               "res_IU_kg_PROPHYLAXIS", "res_IU_kg_joint_bleed",
                               "res_IU_kg_other_bleed", "res_FIX_bleed_ratio",
                               "res_hosp_prob", "res_hosp_LOS",
                               
                               "ETRANACOGENE_abr",
                               "ETRANACOGENE_relative_bleed_reduction",
                               "ETRANACOGENE_max_bleed_reduction_duration",
                               "ETRANACOGENE_bleed_increase_per_year",
                               "ETRANACOGENE_success_prob",
                               
                               "PROPHYLAXIS_abr",
                               "PROPHYLAXIS_relative_bleed_reduction"
                               )
  
  # Create array based on list of parameter names and number of simulations
  
  model$params$prob_sample <- array(0,
                                    dim=c(model$struct$nsim,length(model$params$prob_names)),
                                    dimnames = list(1:model$struct$nsim,
                                                    model$params$prob_names))
  
  # For each simulation, we assign the sampled coefficients to our array of probabilistic parameters
  
  for (i in 1:model$struct$nsim) {
    
    #===============================================================================
    # Population Generation Parameters
    #===============================================================================
    
    #-------------------------------------------------------------------------------
    # Female Share of Patients
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # World Federation of Hemophilia Annual Global Survey (2015)
    # https://www.hemophilia.org/research/research-projects/the-wfh-annual-global-survey-gender-distribution
    
    # DETERMINISTIC:
    # Out of 28'385 reported instances of hemophilia, 1'328 were female
    # That is about 4.7%
    
    # PROBABILISTIC:
    # As we are dealing with a probability (between 0 and 100%) we can sample this directly with a beta distribution
    # We thus do not need to calculate SD and SE 
    # but can directly use the number of females in the sample and the sample size to calculate our distribution
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"share_female"] <- model$params$share_female$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"share_female"] <- rbeta(1,
                                                          shape1 = model$params$share_female$count + 1, 
                                                          shape2 = model$params$share_female$nsample - model$params$share_female$count + 1)
      
    }
    
    #-------------------------------------------------------------------------------
    # Historical Annualized Bleed Rate (ABR)
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Historic real-world ABR of patients.
    # "Real-world comparative analysis of bleeding complications and health-related quality of life in patients with haemophilia A and haemophilia B"
    # Table 2 (Booth et al. 2018)
    
    # DETERMINISTIC:
    # Mean ABR can be implemented directly from the publication.
    
    # PROBABILISTIC:
    # ABR follows a distribution from 0 to positive infinity, which we can model with a gamma distribution.
    # With mean, sample size, and standard deviation, we can calculate the standard error.
    # This can be used to calculate alpha (shape) and beta (scale) parameters of gamme distribution. 
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"hist_abr_exp"] <- model$params$hist_abr$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"hist_abr_exp"] <- rgamma(1, shape = model$params$hist_abr$alpha, scale = model$params$hist_abr$beta)
      
    }
    
    #-------------------------------------------------------------------------------
    # Baseline Annualized Bleed Rate (ABR)
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Baseline annualized bleed rate (ABR) without treatment, based on ABR of patients that received on-demand infusions.
    # Based on assumption that being treated on-demand only treats bleeds but does not lower ABR.
    # 25 on-demand patients were observed for 26 weeks, and had 417 bleeds in this time period. 
    # This can be used to calculate an ABR.
    # "Once-weekly prophylactic treatment vs. on-demand treatment with nonacog alfa in patients with moderately severe to severe haemophilia B" (Kavakli et. al 2016)
    
    # DETERMINISTIC:
    # 25 on-demand patients were observed for 26 weeks, and had 417 bleeds in this time period. 
    # This can be used to calculate an ABR.
    
    # PROBABILISTIC:
    # ABR follows a distribution from 0 to positive infinity, which we can model with a gamma distribution.
    # With mean, sample size, and standard deviation, we can calculate the standard error.
    # This can be used to calculate alpha (shape) and beta (scale) parameters of gamma distribution.
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"baseline_abr_exp"] <- model$params$baseline_abr$mean
        
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"baseline_abr_exp"] <- rgamma(1, shape = model$params$baseline_abr$alpha, scale = model$params$baseline_abr$beta)
      
    }
    
    #-------------------------------------------------------------------------------
    # Historic/Baseline ABR Ratio
    #-------------------------------------------------------------------------------
    
    # Based on the sampled historic and baseline ABR, a ratio is calculated.
    # This guarantees that there is a fixed relationship between historic and baseline ABR for each patient within a simulation.
    # Baseline ABR is sampled individually per patient around the expected probabilistic baseline to account for heterogeneity.
    # If historic and baseline were sampled and applied independently per patient, 
    # we could end up with clinically implausible combinations of historic and baseline ABR within the individual.
    # The Historic/Baseline ABR ratio is then applied to the baseline ABR of every patient.
    
    model$params$prob_sample[i,"hist_baseline_abr_ratio"] <- model$params$prob_sample[i,"hist_abr_exp"]/model$params$prob_sample[i,"baseline_abr_exp"]
    
    #===============================================================================
    # Prices
    #===============================================================================
    
    #-------------------------------------------------------------------------------
    # Etranacogene Dezaparvovec 
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # 
    
    # DETERMINISTIC:
    # 
    
    # PROBABILISTIC:
    # We do not vary ths price probabilistically
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"price_ETRANACOGENE"] <- model$params$prices$ETRANACOGENE$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"price_ETRANACOGENE"] <- model$params$prices$ETRANACOGENE$mean
      
    }
    
    #-------------------------------------------------------------------------------
    # Coagulation factor: Nonacog beta pegol
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Burke et al. (2021) "Clinical, humanistic, and economic burden of severe haemophilia B in adults receiving factor IX prophylaxis: findings from the CHESS II real‑world burden of illness study in Europe"
    
    # DETERMINISTIC:
    # Price Refixia / Nonacog Beta Pegol:  1.70 per IU
    
    # PROBABILISTIC:
    # We do not vary price of the coagulation factor probabilistically
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"price_coagulation_factor_IU"] <- model$params$prices$coagulation_factor_IU$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"price_coagulation_factor_IU"] <- model$params$prices$coagulation_factor_IU$mean
      
    }
    
    #-------------------------------------------------------------------------------
    # Hospitalization: ICU Days
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Ohara et al. (2017): 
    # "The cost of severe haemophilia in Europe: the CHESS study"
    
    # DETERMINISTIC:
    # A day in the ICU after a bleed is estimated to cost EUR 1265.00
    # As this data is from 2014, we adjust for inflation with the CCEMG - EPPI cost converter:
    # https://eppi.ioe.ac.uk/costconversion/default.aspx
    # Cost conversion factor from 2014 to 2022: 1.161164
    
    # PROBABILISTIC:
    # We do not vary the price for a day in the ICU probabilistically
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"price_hosp_ICU_days"] <- model$params$prices$ICU_days$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"price_hosp_ICU_days"] <- model$params$prices$ICU_days$mean
      
    }
    
    #-------------------------------------------------------------------------------
    # Hospitalization: Ward Days
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Ohara et al. (2017): 
    # "The cost of severe haemophilia in Europe: the CHESS study"
    
    # DETERMINISTIC:
    # A day in the hospital ward after a bleed is estimated to cost EUR 514.29
    # As this data is from 2014, we adjust for inflation with the CCEMG - EPPI cost converter:
    # https://eppi.ioe.ac.uk/costconversion/default.aspx
    # Cost conversion factor from 2014 to 2022: 1.161164
    
    # PROBABILISTIC:
    # We do not vary the price for a day in the hospital ward probabilistically
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"price_hosp_ward_days"] <- model$params$prices$hosp_ward_days$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"price_hosp_ward_days"] <- model$params$prices$hosp_ward_days$mean
      
    }
    
    #-------------------------------------------------------------------------------
    # Surgery
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Ohara et al. (2017): 
    # "The cost of severe haemophilia in Europe: the CHESS study"
    
    # DETERMINISTIC:
    # The price of a procedure on a target joint is estimated between 12.02 - 1719.43
    # As we do not know further what determines this range, we take the average of the two and also adjust for inflation
    
    # PROBABILISTIC:
    # We do not vary the price for joint surgery probabilistically
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"price_surgery"] <- model$params$prices$surgery$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"price_surgery"] <- model$params$prices$surgery$mean
      
    }
    
    #===============================================================================
    # Resource Use
    #===============================================================================
    
    #-------------------------------------------------------------------------------
    # Resource Use: IU/Kg of coagulation factor with Prophylaxis (Nonacog Beta Pegol)
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # EMA 2017: https://www.ema.europa.eu/en/documents/product-information/refixia-epar-product-information_en.pdf 
    
    # DETERMINISTIC:
    # 40 IU/kg body weight once weekly
    
    # PROBABILISTIC:
    # No probabilistic component, same as deterministic
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"res_IU_kg_PROPHYLAXIS"] <- model$params$resources$IU_kg_PROPHYLAXIS
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"res_IU_kg_PROPHYLAXIS"] <- model$params$resources$IU_kg_PROPHYLAXIS
      
    }
    
    #-------------------------------------------------------------------------------
    # Resource Use: IU/Kg of Coagulation Factor per joint bleed (Nonacog Beta Pegol)
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # EMA 2017: https://www.ema.europa.eu/en/documents/product-information/refixia-epar-product-information_en.pdf
    
    # DETERMINISTIC:
    # 40 IU/kg body weight once per haemarthrosis
    
    # PROBABILISTIC:
    # No probabilistic component, same as deterministic
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"res_IU_kg_joint_bleed"] <- model$params$resources$IU_kg_joint_bleed
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"res_IU_kg_joint_bleed"] <- model$params$resources$IU_kg_joint_bleed
      
    }
    
    #-------------------------------------------------------------------------------
    # Resource Use: IU/Kg of Coagulation Factor per other bleed (Nonacog Beta Pegol)
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # EMA 2017: https://www.ema.europa.eu/en/documents/product-information/refixia-epar-product-information_en.pdf 
    
    # DETERMINISTIC:
    # 40 IU/kg body weight per non-severe bleed
    
    # PROBABILISTIC:
    # No probabilistic component, same as deterministic
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"res_IU_kg_other_bleed"] <- model$params$resources$IU_kg_other_bleed
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"res_IU_kg_other_bleed"] <- model$params$resources$IU_kg_other_bleed
      
    }
    
    #-------------------------------------------------------------------------------
    # Resource Use: Ratio of bleeds which require FIX infusion
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # "Gene Therapy with Etranacogene Dezaparvovec for Hemophilia B" (Pipe et. al 2023)
    
    # DETERMINISTIC:
    # In total, 148 out of 190 bleeds in the trial required FIX infusion
    
    # PROBABILISTIC:
    # PROBABILISTIC:
    # As we are dealing with a probability (between 0 and 100%) we can sample this directly with a beta distribution
    # We thus do not need to calculate SD and SE 
    # but can directly use the number of bleeds and bleeds requiring FIX infusions to calculate our shape parameters
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"res_FIX_bleed_ratio"] <- model$params$resources$FIX_bleeds_ratio$mean 
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"res_FIX_bleed_ratio"] <- rbeta(1,
                                                                 shape1 = model$params$resources$FIX_bleeds_ratio$FIX_bleeds + 1, 
                                                                 shape2 = model$params$resources$FIX_bleeds_ratio$total_bleeds - model$params$resources$FIX_bleeds_ratio$FIX_bleeds + 1)
    }
    
    #-------------------------------------------------------------------------------
    # Resource Use: Probability of Hospitalization per bleed
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Burke et al. (2021):
    # "Clinical, humanistic, and economic burden
    # of severe haemophilia B in adults receiving factor IX prophylaxis: findings from the CHESS II
    # real‑world burden of illness study in Europe
    
    # DETERMINISTIC:
    # We estimate the ratio of hospitalizations to bleeds: 50 / 180 = 0.27777
    
    # PROBABILISTIC:
    # As we are dealing with a probability (between 0 and 100%) we can sample this directly with a beta distribution
    # We thus do not need to calculate SD and SE 
    # but can directly use the number of bleeds and hospitalizations to calculate our shape parameters
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"res_hosp_prob"] <- model$params$resources$hosp_ratio$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"res_hosp_prob"] <- rbeta(1,
                                                           shape1 = model$params$resources$hosp_ratio$hosps + 1, 
                                                           shape2 = model$params$resources$hosp_ratio$bleeds - model$params$resources$hosp_ratio$hosps + 1)
      
    }
    
    #-------------------------------------------------------------------------------
    # Resource Use: Hospital Length of Stay (LOS)
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Burke et al. (2021):
    # "Clinical, humanistic, and economic burden
    # of severe haemophilia B in adults receiving factor IX prophylaxis: findings from the CHESS II
    # real‑world burden of illness study in Europe
    
    # DETERMINISTIC:
    # Mean bleed-related hospital days per patient: 1.5 (SD = 4.1) based on 62 observed patients
    
    # PROBABILISTIC:
    # Using mean, sample size, and standard deviation, we calculate the standard error
    # We then calculate Alpha and Beta for the Gamma distribution
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"res_hosp_LOS"] <- model$params$resources$hosp_LOS$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"res_hosp_LOS"] <- rgamma(1, shape = model$params$resources$hosp_LOS$alpha, scale = model$params$resources$hosp_LOS$beta)
      
    }
    
    #===============================================================================
    # Utilities
    #===============================================================================
    
    #-------------------------------------------------------------------------------
    # Utilities: Baseline, Age, Sex
    #-------------------------------------------------------------------------------
    
    # Deterministic:
    # We take the estimated coefficients from Ara and Brazier (2010) for all simulations
    
    # Probabilistic:
    # We use the multinorminv function 
    # to randomly draw coefficients based on Cholesky decomposition 
    # of covariance matrix from Ara and Brazier (2010)
    # This allows us to jointly sample coefficients for the baseline utility,
    # as well as the (dis)utility based on sex, age, and age-squared
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"util_male"]     <- model$params$utilities$mu[1]
      model$params$prob_sample[i,"util_age"]      <- model$params$utilities$mu[2]
      model$params$prob_sample[i,"util_agesq"]    <- model$params$utilities$mu[3]
      model$params$prob_sample[i,"util_baseline"] <- model$params$utilities$mu[4]
      
    } else if (model$struct$mode == "probabilistic") {
      
      temp_utils <- fun_multinorminv(model$params$utilities$mu,
                                     model$params$utilities$cov_matrix,
                                     qnorm(runif(4), mean = 0, sd = 1))
      
      model$params$prob_sample[i,"util_male"]     <- temp_utils[1]
      model$params$prob_sample[i,"util_age"]      <- temp_utils[2]
      model$params$prob_sample[i,"util_agesq"]    <- temp_utils[3]
      model$params$prob_sample[i,"util_baseline"] <- temp_utils[4]
      
    }
    
    #-------------------------------------------------------------------------------
    # Utilities: Arthropathy
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Arthropathy disutilities based on "The association of haemophilic arthropathy with Health-Related Quality of Life"
    # Fischer et al. (2016)
    
    # DETERMINISTIC:
    # Mean observed SF-36 utility with Pettersson Score (PS) from 0 - 12 : 0.80 (SD: 0.13)
    # Mean observed SF-36 utility with Pettersson Score (PS) from 13 - 21 : 0.77 (SD: 0.13)
    # Mean observed SF-36 utility with Pettersson Score (PS) that is 22 or higher : 0.70 (SD: 0.12)
    # We calculate the disutilities based on this.
    # The disutilities sum up, so when applying them for a PS higher than 22, we calculate 0.80 - 0.03 - 0.07 = 0.70
    
    # PROBABILISTIC:
    # We use the sum of the variances and the sample sizes to calculate the standard errors.
    # We use the standard errors and sample sizes to calculate alpha (shape) and beta (scale) parameters of the Gamma distribution.
    # We turn the parameter negative as it is a disutility.
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"util_PS_13_to_21"]  <- -model$params$utilities$PS_13_to_21$mean
      model$params$prob_sample[i,"util_PS_22_plus"]   <- -model$params$utilities$PS_22_plus$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"util_PS_13_to_21"]  <- -rgamma(1, shape = model$params$utilities$PS_13_to_21$alpha, scale = model$params$utilities$PS_13_to_21$beta)
      model$params$prob_sample[i,"util_PS_22_plus"]   <- -rgamma(1, shape = model$params$utilities$PS_22_plus$alpha, scale = model$params$utilities$PS_22_plus$beta)
      
    }
    
    #-------------------------------------------------------------------------------
    # Utilities: Bleeding
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Neufeld et al. (2012): # "Effect of Acute Bleeding on Daily Quality of Life Assessments in Patients with 
    # Congenital Hemophilia with Inhibitors and Their Families: Observations from the Dosing Observational Study in Hemophilia"
    
    # DETERMINISTIC:
    # We estimate a utility decrement of 0.2 per day with bleed based on difference between bleeding and non-bleeding days.
    # We turn the parameter negative as it is a disutility.
    # Based on this data, we do not differentiate between types of bleeds in regards to utility, but allow this to be added
    # with further data in the future.
    
    # PROBABILISTIC:
    # We estimate the variance of the utility decrement via the sum of the variance of the two utility states.
    # We then transform the sum of the variance back into the SD, giving us a measure of uncertainty of the utility decrement.
    # The utility decrement follows a distribution from 0 to positive infinity, which we can model with a gamma distribution.
    # With mean, sample size, and standard deviation, we can calculate the standard error.
    # This can be used to calculate alpha (shape) and beta (scale) parameters of the Gamma distribution.
    # As we are sampling for the bleed disutilities separately, despite using the same parameters, we may draw different disutilities from our distribution.
    # We turn the parameter negative as it is a disutility.
    
    if (model$struct$mode == "deterministic") {
      
      # Disutilities for bleed events (Placeholder)
      model$params$prob_sample[i,"util_joint_bleed"] <- -model$params$utilities$bleed$mean
      model$params$prob_sample[i,"util_other_bleed"] <- -model$params$utilities$bleed$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"util_joint_bleed"] <- -rgamma(1, shape = model$params$utilities$bleed$alpha, scale = model$params$utilities$bleed$beta)
      model$params$prob_sample[i,"util_other_bleed"] <- -rgamma(1, shape = model$params$utilities$bleed$alpha, scale = model$params$utilities$bleed$beta)
      
    }
    
    #-------------------------------------------------------------------------------
    # Utilities: Infusion
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # "Preferences and Health-Related Quality-of-Life Related to Disease and Treatment Features for Patients with Hemophilia A in a Canadian General Population Sample" (Johnston et. al 2021)
    
    # DETERMINISTIC:
    # Disutility for infusion with coagulation factor
    # Study shows -0.0003 per infusion per year
    # If prophylaxis includes 52 infusions, 1 per week, this is a total disutility of -0.0003 * 52 = -0.0156 in terms of a health state
    # This disutility is applied to everyone in a health state that includes prophylaxis
    # We turn the parameter negative as it is a disutility.
    
    # PROBABILISTIC:
    # We calculate Alpha and Beta for the Gamma distribution
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"util_infusion"] <- -model$params$utilities$infusion$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"util_infusion"] <- -rgamma(1, shape = model$params$utilities$infusion$alpha, scale = model$params$utilities$infusion$beta)
      
    }
    
    #-------------------------------------------------------------------------------
    # Utilities: Surgery
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Disutility for joint replacement surgery based on Carroll et al. (2019) "Real-world utilities and health-related quality-of life data in hemophilia patients in France and the United Kingdom"
    # Calculation: 0.82 - 0.64 = 0.18
    # Duration of disutility based on Length of Stay (LOS) from Ballal et. al (2008) "Economic evaluation of major knee surgery with recombinant activated factor VII in hemophilia patients with high titer inhibitors and advanced knee arthropathy: exploratory results via literature-based modeling"
    # Length of Stay: 27 Days
    # We scale the disutility based on the length of stay, making the assumption that a cycle is no shorter than 27 days.
    
    # DETERMINISTIC:
    # -0.18 * (27/cycle length)
    
    # PROBABILISTIC:
    # We calculate Alpha and Beta for the Gamma distribution
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"util_surgery"] <- -model$params$utilities$surgery$mean * 27/model$struct$cyc2day
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"util_surgery"] <- -rgamma(1, shape = model$params$utilities$surgery$alpha, scale = model$params$utilities$surgery$beta) * 27/model$struct$cyc2day
      
    }
    
    #===============================================================================
    # Clinical Parameters
    #===============================================================================
    
    #-------------------------------------------------------------------------------
    # Clinical Parameters: Bleeding Type Frequency
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # "Changing patterns of bleeding in patients with severe haemophilia A" (Stephensen et. al 2009) 

    # DETERMINISTIC:
    # Share of joint bleeds: 61.7%
    # Share of other bleeds: 38.3%
    
    # PROBABILISTIC:
    # We use a Beta distribution to draw random proportion of bleeding frequencies.
    # Beta uses the counts of the different categories of bleeds to draw random proportions with a sum of 1.
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"share_joint_bleeds"] <- model$params$joint_bleeds$share
      model$params$prob_sample[i,"share_other_bleeds"] <- model$params$other_bleeds$share
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"share_joint_bleeds"] <- rbeta(1,
                                                                shape1 = model$params$joint_bleeds$events + 1, 
                                                                shape2 =  model$params$joint_bleeds$nsample - model$params$joint_bleeds$events + 1)

      model$params$prob_sample[i,"share_other_bleeds"] <- 1 - model$params$prob_sample[i,"share_joint_bleeds"]
      
    }
    
    #-------------------------------------------------------------------------------
    # Clinical Parameters: Mortality from bleeding
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # "Mortality, life expectancy, and causes of death of persons with hemophilia in the Netherlands 2001–2018" (Hassan et. al 2021) 
    
    # DETERMINISTIC:
    # Standardized Mortality Ratio (SMR) of 2.4 (CI 1.8 - 3.0) for severe hemophilia patients
    # We make the simplifying assumption that this SMR is reached at historic ABR.
    # If ABR is reduced in comparison, SMR also falls.
    # In this manner, more effective treatment can also reduce mortality.
    
    # PROBABILISTIC:
    # Since CI was calculated with normal distribution, we can calculate standard error (SE) from this
    # 95% confidence interval is 3.92 standard errors wide (2 * 1.96)
    # SE = (upper limit - lower limit)/3.92
    # With mean and standard error, we can estimate alpha and beta to fit a Gamma distribution.
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"bleed_mortality_ratio"] <- model$params$bleed_mort$mean

    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"bleed_mortality_ratio"] <- rgamma(1, shape = model$params$bleed_mort$alpha, scale = model$params$bleed_mort$beta)
      
    }
    
    #-------------------------------------------------------------------------------
    # Clinical Parameters: Pettersson Score
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # Association between joint bleeds and Pettersson scores in severe Haemophilia (Fischer et al. 2002)
    
    # DETERMINISTIC:
    # Pettersson Score increased 1 point per 13 joint bleeds
    
    # PROBABILISTIC:
    # We fit a gamma distribution with alpha and beta, estimated from mean and standard error
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"PS_bleeds"] <- model$params$pettersson$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"PS_bleeds"] <- rgamma(1, shape = model$params$pettersson$alpha, scale = model$params$pettersson$beta)
      
    }
    
    #===============================================================================
    # Treatments
    #===============================================================================
    
    #-------------------------------------------------------------------------------
    # Treatment 1: Etranacogene Dezaparvovec
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # "Gene Therapy with Etranacogene Dezaparvovec for Hemophilia B" (Pipe et. al 2023)
  
    # DETERMINISTIC:
    
    # ABR WITH VIRAL GENE THERAPY: Mean ABR of 1.51
    # RELATIVE BLEED REDUCTION: Relative reduction of the ABR from our mean baseline, in percent
    # MAXIMUM BLEED REDUCTION PERIOD: No data, we can only make assumptions
    # BLEED RATE INCREASE PER YEAR: No data, we can only make assumptions
    # TREATMENT SUCCESS PROBABILITY: 52/54 (96.3%) of patients responded to treatment and could discontinue prophylaxis
    
    # PROBABILISTIC:
    
    # ABR WITH VIRAL GENE THERAPY: Gamma distribution
    # RELATIVE BLEED REDUCTION: Not sampled, but calculated same as if deterministic, as relative reduction of the ABR from our mean baseline, in percent
    # MAXIMUM BLEED REDUCTION PERIOD: Uniform distribution from lower to upper limit
    # BLEED RATE INCREASE PER YEAR: Uniform distribution from lower to upper limit
    # TREATMENT SUCCESS PROBABILITY: Beta distribution
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"ETRANACOGENE_abr"]                          <- model$params$ETRANACOGENE$abr$mean
      
      model$params$prob_sample[i,"ETRANACOGENE_relative_bleed_reduction"]     <- model$params$ETRANACOGENE$relative_bleed_reduction$mean
      
      model$params$prob_sample[i,"ETRANACOGENE_max_bleed_reduction_duration"] <- model$params$ETRANACOGENE$max_bleed_reduction_duration$mean
      
      model$params$prob_sample[i,"ETRANACOGENE_bleed_increase_per_year"]      <- model$params$ETRANACOGENE$bleed_increase_per_year$mean 
      
      model$params$prob_sample[i,"ETRANACOGENE_success_prob"]                 <- model$params$ETRANACOGENE$success_prob$mean
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"ETRANACOGENE_abr"]                          <- rgamma(1, 
                                                                                            shape = model$params$ETRANACOGENE$abr$alpha, 
                                                                                            scale = model$params$ETRANACOGENE$abr$beta)
      
      model$params$prob_sample[i,"ETRANACOGENE_relative_bleed_reduction"]     <- 1 - (model$params$prob_sample[i,"ETRANACOGENE_abr"]/model$params$baseline_abr$mean)
      
      model$params$prob_sample[i,"ETRANACOGENE_max_bleed_reduction_duration"] <- runif(1, 
                                                                                           min = model$params$ETRANACOGENE$max_bleed_reduction_duration$l_limit, 
                                                                                           max = model$params$ETRANACOGENE$max_bleed_reduction_duration$u_limit)
      
      model$params$prob_sample[i,"ETRANACOGENE_bleed_increase_per_year"]      <- runif(1, 
                                                                                           min = model$params$ETRANACOGENE$bleed_increase_per_year$l_limit, 
                                                                                           max = model$params$ETRANACOGENE$bleed_increase_per_year$u_limit)
      
      model$params$prob_sample[i,"ETRANACOGENE_success_prob"]                 <- rbeta(1,
                                                                                           shape1 = model$params$ETRANACOGENE$success_prob$events + 1, 
                                                                                           shape2 = model$params$ETRANACOGENE$success_prob$nsample - model$params$ETRANACOGENE$success_prob$events + 1)
    }
    
    #-------------------------------------------------------------------------------
    # Treatment 2: Prophylaxis
    #-------------------------------------------------------------------------------
    
    # SOURCE:
    # "Gene Therapy with Etranacogene Dezaparvovec for Hemophilia B" (Pipe et. al 2023)
   
    # DETERMINISTIC:
    # Annualized bleed rate (ABR) during lead-in period of trial
    # Estimated Rate ABR: 4.19
    # Based on the estimated ABR, we calculate an expected reduction in bleeds from the baseline ABR 
    # and apply this to all patients receiving prophylaxis
    
    # PROBABILISTIC:
    # We use the gamma distribution to sample a probabilistic mean ABR
    # Based on the probabilistic mean ABR, we calculate an expected reduction in bleeds from the probabilistic baseline ABR 
    # and apply this to all patients receiving prophylaxis
    
    if (model$struct$mode == "deterministic") {
      
      model$params$prob_sample[i,"PROPHYLAXIS_abr"]                      <- model$params$PROPHYLAXIS$abr$mean
      
      model$params$prob_sample[i,"PROPHYLAXIS_relative_bleed_reduction"] <- model$params$PROPHYLAXIS$relative_bleed_reduction
      
    } else if (model$struct$mode == "probabilistic") {
      
      model$params$prob_sample[i,"PROPHYLAXIS_abr"]                      <- rgamma(1, shape = model$params$PROPHYLAXIS$abr$alpha, scale = model$params$PROPHYLAXIS$abr$beta)
      
      model$params$prob_sample[i,"PROPHYLAXIS_relative_bleed_reduction"] <-  1 - (model$params$prob_sample[i,"PROPHYLAXIS_abr"]/model$params$baseline_abr$mean)
      
    }
  }
  
  return(model)
  
}

