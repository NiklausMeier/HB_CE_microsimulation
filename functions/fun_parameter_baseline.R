#' Function to set baseline parameters for the model
#' This includes:
#' 1. Fixed parameters
#' 2. Population parameters
#' 3. Probabilistic parameters
#' 4. Treatment parameters

fun_parameter_baseline <- function(model){
  
  ################################################################################
  #                                                                              #
  # Fixed Parameters                                                             #
  #                                                                              #
  ################################################################################
  
  #-------------------------------------------------------------------------------
  # Price: Treatments
  #-------------------------------------------------------------------------------
  
  # We do not know price of Etrenacogene Dezaparvovec will be exactly in Germany.
  # In USA, price tag of approx. 3.5 million USD: https://www.drugdiscoverytrends.com/hemophilia-gene-therapy-hemgenix-sets-record-for-worlds-most-expensive-drug/
  # Price of Roctavian (Valoctocogene roxaparvovec), a gene therapy for haemophilia A also by Biomarin,
  # is expected to cost approximately 1.5 million Euros, net of all discounts, in Europe
  # https://www.biopharmadive.com/news/with-european-approval-secured-biomarin-puts-roughly-15m-price-tag-on-he/630501/
  # We use this as our base case price

  model$params$prices$ETRANACOGENE$mean   <- 1500000
  
  # For the price, we cite Burke et al. (2021) "Clinical, humanistic, and economic burden
  # of severe haemophilia B in adults receiving factor IX prophylaxis: findings from the CHESS II
  # real‑world burden of illness study in Europe"
  # Price Refixia: 40 IU/kg
  
  model$params$prices$coagulation_factor_IU$mean <- 1.70
  
  #-------------------------------------------------------------------------------
  # Price: Hospital and Surgery
  #-------------------------------------------------------------------------------
  
  # Ohara et al. (2017): 
  # "The cost of severe haemophilia in Europe: the CHESS study"
  
  # A day in the hospital ward after a bleed is estimated to cost EUR 514.29
  # A day in the Intensive Care Unit (ICU) after a bleed is estimated to cost EUR 1265.00
  # We assume emergency hospitalizations for joint and other bleeds are treated in the ICU-

  # As this data is from 2014, we adjust for inflation with the CCEMG - EPPI cost converter:
  # https://eppi.ioe.ac.uk/costconversion/default.aspx
  
  # Cost conversion factor from 2014 to 2022: 1.161164
  
  model$params$cost_conversion_2014_2022     <- 1.161164
  
  model$params$prices$hosp_ward_days$mean    <- 514.29 * model$params$cost_conversion_2014_2022
  model$params$prices$ICU_days$mean          <- 1265.00 * model$params$cost_conversion_2014_2022

  # The price of a procedure on a target joint is estimated between 12.02 - 1719.43
  # As we do not know further what determines this range, we take the average of the two and also adjust for inflation
  
  model$params$prices$surgery$mean           <- (12.02+1719.43)/2 * model$params$cost_conversion_2014_2022
  
  #-------------------------------------------------------------------------------
  # Background Mortality
  #-------------------------------------------------------------------------------
  
  ## Calculate death rates and probabilities for a given age based on chosen model and cycle length:
  
  model$life_tables_weight <- life_tables_weight
  
  # Calculate death rate for the chosen cycle length during a particular age:
  
  for (i in 0:nrow(life_tables_weight)){
    
    # Female
    model[["life_tables_weight"]][i,"female_cycle_death_rate"]        <- model[["life_tables_weight"]][["female_yearly_death_rate"]][i] * (model$struct$cyc2day/settings$time$yr2day)
    model[["life_tables_weight"]][i,"female_cycle_death_probability"] <- 1-exp(-model[["life_tables_weight"]][["female_cycle_death_rate"]][i])
    
    # Male
    model[["life_tables_weight"]][i,"male_cycle_death_rate"]        <- model[["life_tables_weight"]][["male_yearly_death_rate"]][i] * (model$struct$cyc2day/settings$time$yr2day)
    model[["life_tables_weight"]][i,"male_cycle_death_probability"] <- 1-exp(-model[["life_tables_weight"]][["male_cycle_death_rate"]][i])
    
  }
  
  #-------------------------------------------------------------------------------
  # Discounting
  #-------------------------------------------------------------------------------
  
  # Set discount rate
  model$params$disc_rate    <- 0.03
  
  # Calculate discount factors
  
  model$disc_factors                          <- data.table(matrix(as.double(NA), nrow = model$struct$ncycle , ncol = 2, 
                                                                   dimnames = list(0:(model$struct$ncycle-1), c("years","disc_factors"))))
  
  for (i in 1:(model$struct$ncycle)){
    
    # Calculate discount factor for a particular cycle in the model
    model[["disc_factors"]][i, "years"        := ((i-1)*(model$struct$cyc2day/settings$time$yr2day)) ]
    model[["disc_factors"]][i, "disc_factors" := (1/(1 + model$params$disc_rate))^((i-1)*(model$struct$cyc2day/settings$time$yr2day)) ]
  }
  
  #-------------------------------------------------------------------------------
  # Resource Use: Coagulation Factor
  #-------------------------------------------------------------------------------
  
  # Dosage of coagulation factor Nonacog Beta Pegol based on EMA 2017:
  # https://www.ema.europa.eu/en/documents/product-information/refixia-epar-product-information_en.pdf 
  
  # Prophylaxis: 40 IU/kg body weight once weekly
  model$params$resources$IU_kg_PROPHYLAXIS <- 40
  
  # 40 IU/kg body weight per non-severe bleed
  model$params$resources$IU_kg_joint_bleed <- 40

  # 40 IU/kg body weight per non-severe bleed  
  model$params$resources$IU_kg_other_bleed <- 40
  
  ################################################################################
  #                                                                              #
  # Population Parameters                                                        #
  #                                                                              #
  ################################################################################
  
  #-------------------------------------------------------------------------------
  # Sex
  #-------------------------------------------------------------------------------
  
  # Share of female population
  # https://www.hemophilia.org/research/research-projects/the-wfh-annual-global-survey-gender-distribution
  # Out of 28'385 reported instances of hemophilia, 1'328 were female
  # That is about 4.7%
  
  model$params$share_female$count       <- 1328
  model$params$share_female$nsample     <- 28385
  
  model$params$share_female$mean       <- model$params$share_female$count / model$params$share_female$nsample
    
  #-------------------------------------------------------------------------------
  # Age
  #-------------------------------------------------------------------------------
  
  # Mean and standard deviation stem from: 
  # "Real-world comparative analysis of bleeding complications and health-related quality of life in patients with haemophilia A and haemophilia B"
  # Table 1 (Booth et al. 2018)
  model$params$age$mean         <- 36.3
  model$params$age$sd           <- 15.3
  
  #-------------------------------------------------------------------------------
  # Historical Annualized Bleed Rate (ABR)
  #-------------------------------------------------------------------------------
  
  # Historic annualized bleed rate (ABR) prior to model, assuming various prophylactic treatments thus far
  # Mean, standard deviation and sample size stem from: 
  # "Real-world comparative analysis of bleeding complications and health-related quality of life in patients with haemophilia A and haemophilia B"
  # Table 2 (Booth et al. 2018)
  # Standard error can be calculated from standard deviation and sample size (278 patients)
  
  model$params$hist_abr$mean         <- 4.6
  model$params$hist_abr$sd           <- 5.8
  model$params$hist_abr$nsample      <- 278
  model$params$hist_abr$se           <- model$params$hist_abr$sd/sqrt(model$params$hist_abr$nsample)
  
  model$params$hist_abr$alpha        <- (model$params$hist_abr$mean/model$params$hist_abr$se)^2
  model$params$hist_abr$beta         <- (model$params$hist_abr$se^2)/model$params$hist_abr$mean
  
  #-------------------------------------------------------------------------------
  # Baseline Annualized Bleed Rate (ABR)
  #-------------------------------------------------------------------------------
  
  # Baseline annualized bleed rate (ABR) without treatment, based on ABR of patients that received on-demand infusions.
  # Based on assumption that being treated on-demand only treats bleeds but does not lower ABR.
  # 25 on-demand patients were observed for 26 weeks, and had 417 bleeds in this time period. 
  # This can be used to calculate an ABR.
  # "Once-weekly prophylactic treatment vs. on-demand treatment with nonacog alfa in patients with moderately severe to severe haemophilia B" (Kavakli et. al 2016)
  
  model$params$baseline_abr$mean         <- 32.9
  model$params$baseline_abr$sd           <- 17.4
  model$params$baseline_abr$nsample      <- 25
  model$params$baseline_abr$se           <- model$params$baseline_abr$sd/sqrt(model$params$baseline_abr$nsample)
  
  model$params$baseline_abr$alpha        <- (model$params$baseline_abr$mean/model$params$baseline_abr$se)^2
  model$params$baseline_abr$beta         <- (model$params$baseline_abr$se^2)/model$params$baseline_abr$mean
  
  ################################################################################
  #                                                                              #
  #  Probabilistic Parameters                                                    #
  #                                                                              #
  ################################################################################
  
  #-------------------------------------------------------------------------------
  # Bleeding type frequency: 
  #-------------------------------------------------------------------------------
  
  # We use a Beta distribution to draw random proportion of bleeding frequencies. 
  # This means we only need counts of the different bleeding types: Joint bleeds and other bleeds
  
  # Source: Stephensen et. al (2009) "Changing patterns of bleeding in patients with severe haemophilia A"
  # Reports share of joint bleeds as being 61.7%
  # We turn this back into count data for our parameter estimation
  # According to study there were 8506 treatments, so we calculate that 61.7% of these treatments were joint bleeds
  
  # Joint bleeds:
  
  model$params$joint_bleeds$events   <- 8506*0.617
  model$params$joint_bleeds$nsample  <- 8506
  model$params$joint_bleeds$share    <- model$params$joint_bleeds$events / model$params$joint_bleeds$nsample
  
  ## Share of other bleeds is the remaining bleeds that are not ICH or joint
  
  model$params$other_bleeds$share   <- 1 - model$params$joint_bleeds$share
  
  #-------------------------------------------------------------------------------
  # Mortality from bleeding
  #-------------------------------------------------------------------------------
  
  # Age-adjusted standardized mortality higher for hemophilia patients compared to general population
  # Source: Hassan et. al (2021) "Mortality, life expectancy, and causes of death of persons with hemophilia in the Netherlands 2001–2018"
  # Standardized Mortality Ratio (SMR) of 2.4 (CI 1.8 - 3.0) for severe hemophilia patients
  # We make the simplifying assumption that this SMR is reached at historic ABR.
  # If ABR is reduced in comparison, SMR also falls.
  # In this manner, more effective treatment can also reduce mortality.
  
  # Since CI was calculated with normal distribution, we can calculate standard error (SE) from this
  # 95% confidence interval is 3.92 standard errors wide (2 * 1.96)
  # SE = (upper limit - lower limit)/3.92
  
  model$params$bleed_mort$mean    <- 2.4   # Mortality is increased by 2.4 at historic ABR
  model$params$bleed_mort$u_limit <- 3.0   # Upper limit
  model$params$bleed_mort$l_limit <- 1.8   # Lower limit
  model$params$bleed_mort$se      <- (model$params$bleed_mort$u_limit - model$params$bleed_mort$l_limit)/3.92
  
  # With mean and standard error, we can estimate alpha and beta to fit a Gamma distribution.
  
  model$params$bleed_mort$alpha  <- (model$params$bleed_mort$mean/model$params$bleed_mort$se)^2
  model$params$bleed_mort$beta   <- (model$params$bleed_mort$se^2)/model$params$bleed_mort$mean
  
  #-------------------------------------------------------------------------------
  # Arthropathy
  #-------------------------------------------------------------------------------
  
  # ASSOCIATION BETWEEN JOINT BLEEDS AND PETTERSSON SCORES IN SEVERE HAEMOPHILIA (Fischer et al. 2002)
  # Maximum of Pettersson Score is 78
  # According to Fischer et al., Pettersson Score increased 1 point per 13 joint bleeds (95% CI 11 - 15)
  # Assuming CI was calculated with normal distribution, we can calculate standard error (SE) from this
  # 95% confidence interval is 3.92 standard errors wide (2 * 1.96)
  # SE = (upper limit - lower limit)/3.92
  
  model$params$pettersson$mean    <- 13   # Pettersson score increases by 1 per 13 bleeds
  model$params$pettersson$u_limit <- 15   # Upper limit
  model$params$pettersson$l_limit <- 11   # Lower limit
  model$params$pettersson$se      <- (model$params$pettersson$u_limit - model$params$pettersson$l_limit)/3.92
  
  # Maximum Pettersson score is 78
  model$params$pettersson$max    <- 78
  
  # With mean and standard error, we can estimate alpha and beta to fit a gamma distribution
  
  model$params$pettersson$alpha  <- (model$params$pettersson$mean/model$params$pettersson$se)^2
  model$params$pettersson$beta   <- (model$params$pettersson$se^2)/model$params$pettersson$mean
  
  # Clinically relevant threshold for Pettersson score is 28  (Fischer et al. 2011)
  # "A modeling approach to evaluate long-term outcome of prophylactic
  # and on demand treatment strategies for severe hemophilia A"
  
  model$params$pettersson$clin_rel    <- 28
  
  #-------------------------------------------------------------------------------
  # Resource Use: Share of bleeds requiring FIX infusion
  #-------------------------------------------------------------------------------
  
  # "Gene Therapy with Etranacogene Dezaparvovec for Hemophilia B" (Pipe et. al 2023)
  
  # In total, 148 out of 190 bleeds in the trial required FIX infusion
  
  model$params$resources$FIX_bleeds_ratio$total_bleeds <- 190
  model$params$resources$FIX_bleeds_ratio$FIX_bleeds   <- 148
  
  model$params$resources$FIX_bleeds_ratio$mean          <- model$params$resources$FIX_bleeds_ratio$FIX_bleeds/model$params$resources$FIX_bleeds_ratio$total_bleeds
  
  #-------------------------------------------------------------------------------
  # Resource Use: Probability of Hospitalization per bleed
  #-------------------------------------------------------------------------------
  
  # We estimate the frequency of hospitalizations per bleed based on Burke et al. (2021):
  # "Clinical, humanistic, and economic burden
  # of severe haemophilia B in adults receiving factor IX prophylaxis: findings from the CHESS II
  # real‑world burden of illness study in Europe"
  
  # Sample size = 75 patients
  # Mean ABR: 2.4 => 75 * 2.4 = 180 observed bleeds
  # Mean number of bleed-related hospitalizations: 0.67 => 75 * 0.67 = 50 bleed-related hospitalizations
  # Ratio: 50/180 = 0.27777  probability of hospitalization per bleed 
  
  model$params$resources$hosp_ratio$bleeds <- 180
  model$params$resources$hosp_ratio$hosps  <- 50
  
  model$params$resources$hosp_ratio$mean   <- model$params$resources$hosp_ratio$hosps/model$params$resources$hosp_ratio$bleeds
  
  # As we are dealing with a probability (between 0 and 100%) we can sample this directly with a beta distribution
  # We thus do not need to calculate SD and SE
  
  #-------------------------------------------------------------------------------
  # Resource Use: Hospital Length of Stay (LOS)
  #-------------------------------------------------------------------------------
  
  # We estimate the bleed-related hospital days per patient based on Burke et al. (2021):
  # "Clinical, humanistic, and economic burden
  # of severe haemophilia B in adults receiving factor IX prophylaxis: findings from the CHESS II
  # real‑world burden of illness study in Europe"
  
  # Mean bleed-related hospital days per patient: 1.5 (SD = 4.1) based on 62 observed patients
  
  model$params$resources$hosp_LOS$mean     <- 1.5
  model$params$resources$hosp_LOS$sd       <- 4.1
  model$params$resources$hosp_LOS$nsample  <- 62
  
  model$params$resources$hosp_LOS$se       <- model$params$resources$hosp_LOS$sd/sqrt(model$params$resources$hosp_LOS$nsample)
  
  # We calculate Alpha and Beta for the Gamma distribution
  
  model$params$resources$hosp_LOS$alpha    <- (model$params$resources$hosp_LOS$mean/model$params$resources$hosp_LOS$se)^2
  model$params$resources$hosp_LOS$beta     <- (model$params$resources$hosp_LOS$se^2)/model$params$resources$hosp_LOS$mean
  
  #-------------------------------------------------------------------------------
  # Resource Use: Surgery
  #-------------------------------------------------------------------------------
  
  # Length of Stay (LOS) from Ballal et. al (2008) "Economic evaluation of major knee surgery with recombinant activated factor VII in hemophilia patients with high titer inhibitors and advanced knee arthropathy: exploratory results via literature-based modeling"
  # Length of Stay: 27 Days
  # When calculating the costs, we apply the costs per hospital day based on the length of stay
  
  model$params$resources$surg_LOS <- 27
  
  #===============================================================================
  # Utilities
  #===============================================================================
  
  #-------------------------------------------------------------------------------
  # BASELINE, AGE- AND SEX-BASED UTILITIES
  #-------------------------------------------------------------------------------
  
  # Ara and Brazier (2010):
  # "Populating an Economic Model with Health State Utility Values: Moving toward Better Practice"
  model$params$utilities$baseline   <- 0.9508566
  model$params$utilities$male       <- 0.0212126
  model$params$utilities$age_coef   <- -0.0002587
  model$params$utilities$agesq_coef <- -0.0000332
  
  # Coefficients and covariance matrix for multivariate sampling in probabilistic analysis
  model$params$utilities$mu         <-       c(0.0212126,-0.0002587,-0.0000332,0.9508566)
  
  model$params$utilities$cov_matrix <- rbind(c( 0.0000071,      0,            0,              0),
                                             c(-0.000000038,    0.00000014,   0,              0),
                                             c( 0.00000000038, -0.0000000015, 0.000000000016, 0),
                                             c(-0.0000025,     -0.0000028,    0.000000028,    0.000061)) 
  
  #-------------------------------------------------------------------------------
  # SCENARIO GROCHTDREIS: BASELINE, AGE- AND SEX-BASED UTILITIES
  #-------------------------------------------------------------------------------
  
  # Grochtdreis et al. (2019):
  # For a scenario, we check if HRQoL changes if we use population means from a 
  # German EQ-5D-5L value set instead.
  # We turn the mean values, based on sex and age bracket, into a baseline
  # utility and a disutility
  
  model$params$scen_utilities$male_baseline      <- 0.94
  model$params$scen_utilities$male_disutil_25_34 <- 0.01
  model$params$scen_utilities$male_disutil_35_44 <- 0.01
  model$params$scen_utilities$male_disutil_45_54 <- 0.03
  model$params$scen_utilities$male_disutil_55_64 <- 0.02
  model$params$scen_utilities$male_disutil_65_74 <- 0.01
  model$params$scen_utilities$male_disutil_75    <- 0.02
  
  model$params$scen_utilities$male_baseline - model$params$scen_utilities$male_disutil_25_34 - model$params$scen_utilities$male_disutil_35_44 - model$params$scen_utilities$male_disutil_45_54 - model$params$scen_utilities$male_disutil_55_64 - model$params$scen_utilities$male_disutil_65_74 - model$params$scen_utilities$male_disutil_75
  
  model$params$scen_utilities$female_baseline      <- 0.94
  model$params$scen_utilities$female_disutil_25_34 <- 0.02
  model$params$scen_utilities$female_disutil_35_44 <- 0.04
  model$params$scen_utilities$female_disutil_45_54 <- 0.02
  model$params$scen_utilities$female_disutil_55_64 <- 0.00
  model$params$scen_utilities$female_disutil_65_74 <- 0.01
  model$params$scen_utilities$female_disutil_75    <- 0.08
  
  model$params$scen_utilities$female_baseline - model$params$scen_utilities$female_disutil_25_34 - model$params$scen_utilities$female_disutil_35_44 - model$params$scen_utilities$female_disutil_45_54 - model$params$scen_utilities$female_disutil_55_64 - model$params$scen_utilities$female_disutil_65_74 - model$params$scen_utilities$female_disutil_75
  
  #-------------------------------------------------------------------------------
  # DISUTILITY ARTHROPATHY
  #-------------------------------------------------------------------------------
  
  # Fischer et al. (2016):
  # "The association of haemophilic arthropathy with Health-Related Quality of Life"
  # Mean observed SF-36 utility with Pettersson Score (PS) from 0 - 12 : 0.80 (SD: 0.13)
  # Mean observed SF-36 utility with Pettersson Score (PS) from 13 - 21 : 0.77 (SD: 0.13)
  # Mean observed SF-36 utility with Pettersson Score (PS) that is 22 or higher : 0.70 (SD: 0.12)
  # We calculate the disutilities based on this
  # We use positive values for the parameters to allow for sampling via the Gamma distribution
  
  model$params$utilities$PS_13_to_21$mean    <- 0.03
  model$params$utilities$PS_22_plus$mean     <- 0.07
  
  # We transform the standard deviations into variances
  model$params$utilities$PS_0_to_12$sd       <- 0.13
  model$params$utilities$PS_13_to_21$sd      <- 0.13
  model$params$utilities$PS_22_plus$sd       <- 0.12
  
  model$params$utilities$PS_0_to_12$var      <- model$params$utilities$PS_0_to_12$sd^2
  model$params$utilities$PS_13_to_21$var     <- model$params$utilities$PS_13_to_21$sd^2
  model$params$utilities$PS_22_plus$var      <- model$params$utilities$PS_22_plus$sd^2
  
  # We need sample sizes to calculate the standard error (SE)
  model$params$utilities$PS_0_to_12$nsample  <- 68
  model$params$utilities$PS_13_to_21$nsample <- 36
  model$params$utilities$PS_22_plus$nsample  <- 72
  
  # With the sum of the variances and the sample sizes (N), 
  # we can calculate the standard error of the difference (SE_diff) with the following formula: 
  # SE_diff = sqrt(var1/N1 + var2/N2)
  
  model$params$utilities$PS_13_to_21$se      <- sqrt(model$params$utilities$PS_0_to_12$var/(model$params$utilities$PS_0_to_12$nsample) + model$params$utilities$PS_13_to_21$var/(model$params$utilities$PS_13_to_21$nsample))
  model$params$utilities$PS_22_plus$se       <- sqrt(model$params$utilities$PS_13_to_21$var/(model$params$utilities$PS_13_to_21$nsample) + model$params$utilities$PS_22_plus$var/(model$params$utilities$PS_22_plus$nsample))
  
  # We calculate Alpha and Beta for the Gamma distribution
  
  model$params$utilities$PS_13_to_21$alpha   <- (model$params$utilities$PS_13_to_21$mean/model$params$utilities$PS_13_to_21$se)^2
  model$params$utilities$PS_13_to_21$beta    <- (model$params$utilities$PS_13_to_21$se^2)/ model$params$utilities$PS_13_to_21$mean
  
  model$params$utilities$PS_22_plus$alpha    <- (model$params$utilities$PS_22_plus$mean/model$params$utilities$PS_22_plus$se)^2
  model$params$utilities$PS_22_plus$beta     <- (model$params$utilities$PS_22_plus$se^2)/ model$params$utilities$PS_22_plus$mean
  
  #-------------------------------------------------------------------------------
  # DISUTILITY BLEEDS
  #-------------------------------------------------------------------------------
  
  # Neufeld et al. (2012):
  # "Effect of Acute Bleeding on Daily Quality of Life Assessments in
  # Patients with Congenital Hemophilia with Inhibitors and Their
  # Families: Observations from the Dosing Observational Study in Hemophilia"
  # Observed EQ-5D utility on days without bleed: 0.84 (SD: 0.16)
  # Observed EQ-5D utility on days with bleed: 0.64 (SD: 0.23)
  # Estimated utility decrement: 0.64 - 0.84 = -0.20 (for one day) based on a sample size of 37 patients
  # We use a positive value for the parameter to allow for sampling via the Gamma distribution
  model$params$utilities$bleed$mean    <- 0.20
  
  # We transform the standard deviations into variances
  model$params$utilities$bleed$sd1     <- 0.16
  model$params$utilities$bleed$var1    <- model$params$utilities$bleed$sd1^2
  
  model$params$utilities$bleed$sd2     <- 0.23
  model$params$utilities$bleed$var2    <- model$params$utilities$bleed$sd2^2
  
  # We estimate the variance of the utility decrement via the sum of the variance of the two utility states
  # With this sum of the variance and the sample sizes (N), we can calculate the standard error of the difference (SE_diff) with the following formula: 
  # SE_diff = sqrt(var1/N1 + var2/N2)
  
  model$params$utilities$bleed$nsample <- 37
  model$params$utilities$bleed$se      <- sqrt(model$params$utilities$bleed$var1/(model$params$utilities$bleed$nsample) + model$params$utilities$bleed$var2/(model$params$utilities$bleed$nsample))
  
  # We calculate Alpha and Beta for the Gamma distribution
  
  model$params$utilities$bleed$alpha   <- (model$params$utilities$bleed$mean/model$params$utilities$bleed$se)^2
  model$params$utilities$bleed$beta    <- (model$params$utilities$bleed$se^2)/model$params$utilities$bleed$mean
  
  #-------------------------------------------------------------------------------
  # DISUTILITY INFUSION
  #-------------------------------------------------------------------------------
  
  # "Preferences and Health-Related Quality-of-Life Related to Disease and Treatment Features for Patients with Hemophilia A in a Canadian General Population Sample" (Johnston et. al 2021)
  # Study shows -0.0003 per infusion per year
  # If prophylaxis includes 52 infusions, 1 per week, this is a total disutility of -0.0003 * 52 = -0.0156 in terms of a health state
  # This disutility is applied to everyone in a health state that includes prophylaxis
  
  model$params$utilities$infusion$mean    <- 0.0156
  
  # We use the p-value (no more than 0.01) to estimate the z-value, with which the standard error (se) can also be calculated
  
  model$params$utilities$infusion$p       <- 0.001
  model$params$utilities$infusion$z       <- (-0.862 + sqrt(0.743 - 2.404*log(model$params$utilities$infusion$p)))
  model$params$utilities$infusion$se      <- (model$params$utilities$infusion$mean)/model$params$utilities$infusion$z
  
  # We calculate Alpha and Beta for the Gamma distribution
  
  model$params$utilities$infusion$alpha   <- (-model$params$utilities$infusion$mean/model$params$utilities$infusion$se)^2
  model$params$utilities$infusion$beta    <- (model$params$utilities$infusion$se^2)/(model$params$utilities$infusion$mean)
  
  #-------------------------------------------------------------------------------
  # DISUTILITY JOINT REPLACEMENT SURGERY
  #-------------------------------------------------------------------------------
  
  # Disutility for joint replacement surgery based on Carroll et al. (2019) "Real-world utilities and health-related quality-of life data in hemophilia patients in France and the United Kingdom"
  # Mean EQ-5D utility with surgery: 0.64 (SD: 0.25)
  # Mean EQ-5D utility without surgery: 0.82 (SD: 0.22)
  # Calculation: 0.82 - 0.64 = 0.18
  # Duration of disutility based on Length of Stay (LOS) from Ballal et. al (2008) "Economic evaluation of major knee surgery with recombinant activated factor VII in hemophilia patients with high titer inhibitors and advanced knee arthropathy: exploratory results via literature-based modeling"
  # Length of Stay: 27 Days
  # We scale the disutility based on the length of stay, making the assumption that a cycle is no shorter than 27 days.
  
  # We use a positive value for the parameter to allow for sampling via the Gamma distribution
  model$params$utilities$surgery$mean        <- 0.18
  
  # We transform the standard deviations into variances
  
  model$params$utilities$surgery$sd1  <- 0.25
  model$params$utilities$surgery$var1 <- model$params$utilities$surgery$sd1^2
    
  model$params$utilities$surgery$sd2  <- 0.22
  model$params$utilities$surgery$var2 <- model$params$utilities$surgery$sd2^2
  
  # Sample Size
  # We have:
  # 1. Sample size no surgery: 60
  # 2. Sample size surgery: 70
  # This gives us a combined sample size of 130
  
  model$params$utilities$surgery$nsample <- 130
  
  # We estimate the variance of the utility decrement via the sum of the variance of the two utility states
  # With this sum of the variance and the sample sizes (N), we can calculate the standard error of the difference (SE_diff) with the following formula: 
  # SE_diff = sqrt(var1/N1 + var2/N2)
  
  model$params$utilities$surgery$se      <- sqrt(model$params$utilities$surgery$var1/(model$params$utilities$surgery$nsample) + model$params$utilities$surgery$var2/(model$params$utilities$surgery$nsample))
  
  # We calculate Alpha and Beta for the Gamma distribution
  
  model$params$utilities$surgery$alpha   <- (-model$params$utilities$surgery$mean/model$params$utilities$surgery$se)^2
  model$params$utilities$surgery$beta    <- (model$params$utilities$surgery$se^2)/(model$params$utilities$surgery$mean)
  
  ################################################################################
  #                                                                              #
  #  TREATMENTS                                                                  #
  #                                                                              #
  ################################################################################
  
  #' We number the treatments so we can track them in our array.
  #' We have three categories of treatment based on the treatment arm:
  #' Treatments 1X: Etranacogene Dezaparvovec
  #' Treatments 2: Prophylaxis

  # More specifically:
  
  #' Treatment 11: Etranacogene Dezaparvovec treatment period success
  #' Treatment 12: Etranacogene Dezaparvovec treatment period failure
  #' Treatment 13: Etranacogene Dezaparvovec observation period
  #' Treatment 14: Etranacogene Dezaparvovec treatment with prophylaxis

  #' Treatment 21: Prophylaxis treatment
  
  #-------------------------------------------------------------------------------
  # Treatment 1: Etranacogene Dezaparvovec
  #-------------------------------------------------------------------------------
  
  # Annualized bleed rate (ABR) while being treated with Etrenacogene Dezaparvovec and being observed for 52 weeks:
  # "Gene Therapy with Etranacogene Dezaparvovec for Hemophilia B" (Pipe et. al 2023)
  
  ## ABR WITH ETRANACOGENE
  # Based on report, we know that 20/54 patients reported 54 bleeds in total during the observation period
  # However, different individuals had different observation periods, so the ABR cannot be calculated directly from this
  # An ABR of 1.51 is reported, with a 95%-CI from a negative binomial regression model of 0.81-2.82

  model$params$ETRANACOGENE$abr$mean     <- 1.51

  # We know that CI was not calculated with normal distribution, but based on a negative binomial regression
  # However, we do not have an SD, and can not recreate the CI from the negative binomial regression
  # Therefore, we act as if CI was calculated with normal distribution, so we can calculate standard error (SE) from this
  # This is not completely accurate, but still helps to capture uncertainty
  # 95% confidence interval is 3.92 standard errors wide (2 * 1.96)
  # SE = (upper limit - lower limit)/3.92
  
  # The negative binomial regression CI is 0.81 - 2.82
  # However, to ensure our mean value is the same in the probabilistic as in the deterministic model, 
  # we must adjust these slightly so that they average out to 1.51
  # (2.82 + 0.81) / 2 = 3.63 / 2 = 1.815
  # 1.815 - 1.51 = 0.305
  # We lower both lower and upper bound of the CI by 0.305
  
  model$params$ETRANACOGENE$abr$u_limit         <- 2.82 - 0.305
  model$params$ETRANACOGENE$abr$l_limit         <- 0.81 - 0.305
  
  model$params$ETRANACOGENE$abr$se              <- (model$params$ETRANACOGENE$abr$u_limit - model$params$ETRANACOGENE$abr$l_limit)/3.92
  
  # Alpha and Beta parameters for Gamma distribution
  
  model$params$ETRANACOGENE$abr$alpha   <- (model$params$ETRANACOGENE$abr$mean/model$params$ETRANACOGENE$abr$se)^2
  model$params$ETRANACOGENE$abr$beta    <- (model$params$ETRANACOGENE$abr$se^2)/model$params$ETRANACOGENE$abr$mean
  
  # Test whether this replicates CI well
  mean(rgamma(10000, shape = model$params$ETRANACOGENE$abr$alpha, scale = model$params$ETRANACOGENE$abr$beta))
  quantile(rgamma(10000, shape = model$params$ETRANACOGENE$abr$alpha, scale = model$params$ETRANACOGENE$abr$beta), c(.025, .500, 0.975))
  
  # The CI is overall a bit lower than the original 
  
  # The median is close to the mean, but the The CI is overall a bit lower than the original 
  # Given the non-symmetric nature of the original CI, some level of inaccuracy is inevitable.
  # The mean and median are roughly appropriate.
  
  # RELATIVE BLEED REDUCTION:
  # We can also express this as a relative reduction of the ABR from our mean baseline, in percent
  
  model$params$ETRANACOGENE$relative_bleed_reduction$mean   <- 1 - (model$params$ETRANACOGENE$abr$mean/model$params$baseline_abr$mean)
  
  ## MAXIMUM BLEED REDUCTION PERIOD
  # We have no data on how long this maximum reduction on ABR will hold 
  # We can only make assumptions
  # We specify an expected value, an upper limit and a lower limit
  
  model$params$ETRANACOGENE$max_bleed_reduction_duration$mean    <- 10
  model$params$ETRANACOGENE$max_bleed_reduction_duration$l_limit <- 5
  model$params$ETRANACOGENE$max_bleed_reduction_duration$u_limit <- 15
  
  ## BLEED RATE INCREASE PER YEAR
  # We have no data on how quickly the bleed rate will again increase over time
  # We can only make assumptions
  # We specify an expected value, an upper limit and a lower limit
  model$params$ETRANACOGENE$bleed_increase_per_year$mean    <- 0.10
  model$params$ETRANACOGENE$bleed_increase_per_year$l_limit <- 0.05
  model$params$ETRANACOGENE$bleed_increase_per_year$u_limit <- 0.15
  
  ## TREATMENT SUCCESS PROBABILITY
  # 52/54 (96.3%) of patients responded to treatment and could discontinue prophylaxis
  # We can fit a beta distribution directly based on the number of events (52) 
  # out of the total sample minus the events (54 - 52 = 2)
  
  model$params$ETRANACOGENE$success_prob$events  <- 52
  model$params$ETRANACOGENE$success_prob$nsample <- 54
  model$params$ETRANACOGENE$success_prob$mean    <- model$params$ETRANACOGENE$success_prob$events/model$params$ETRANACOGENE$success_prob$nsample
  
  #-------------------------------------------------------------------------------
  # Treatment 2: Prophylaxis
  #-------------------------------------------------------------------------------
  
  # Annualized bleed rate (ABR) during lead-in period before treatment with Etranacogene Dezaparvovec
  # "Gene Therapy with Etranacogene Dezaparvovec for Hemophilia B" (Pipe et. al 2023)
  
  ## ABR during lead-in period: 4.19
  
  model$params$PROPHYLAXIS$abr$mean             <- 4.19
  
  # We know that CI was not calculated with normal distribution, but based on a negative binomial regression
  # However, we do not have an SD, and can not recreate the CI from the negative binomial regression
  # Therefore, we act as if CI was calculated with normal distribution, so we can calculate standard error (SE) from this
  # This is not completely accurate, but still helps to capture uncertainty
  # 95% confidence interval is 3.92 standard errors wide (2 * 1.96)
  # SE = (upper limit - lower limit)/3.92
  
  # The negative binomial regression CI is 3.22 - 5.45
  # However, to ensure our mean value is the same in the probabilistic as in the deterministic model, 
  # we must adjust these slightly so that they average out to 4.19
  # (5.45 + 3.22) / 2 = 8.67 / 2 = 4.335
  # 4.335 - 4.19 = 0.145
  # We lower both lower and upper bound of the CI by 0.145
  
  model$params$PROPHYLAXIS$abr$u_limit         <- 5.45 - 0.145
  model$params$PROPHYLAXIS$abr$l_limit         <- 3.22 - 0.145
  
  model$params$PROPHYLAXIS$abr$se              <- (model$params$PROPHYLAXIS$abr$u_limit - model$params$PROPHYLAXIS$abr$l_limit)/3.92
  
  # Alpha and Beta parameters for Gamma distribution
  
  model$params$PROPHYLAXIS$abr$alpha           <- (model$params$PROPHYLAXIS$abr$mean/model$params$PROPHYLAXIS$abr$se)^2
  model$params$PROPHYLAXIS$abr$beta            <- (model$params$PROPHYLAXIS$abr$se^2)/model$params$PROPHYLAXIS$abr$mean
  
  # RELATIVE BLEED REDUCTION:
  # We can also express this as a relative reduction of the ABR from our mean baseline, in percent
  
  model$params$PROPHYLAXIS$relative_bleed_reduction   <- 1 - (model$params$PROPHYLAXIS$abr$mean/model$params$baseline_abr$mean)
  
  # Test whether this replicates CI well
  mean(rgamma(10000, shape = model$params$PROPHYLAXIS$abr$alpha, scale = model$params$PROPHYLAXIS$abr$beta))
  quantile(rgamma(10000, shape = model$params$PROPHYLAXIS$abr$alpha, scale = model$params$PROPHYLAXIS$abr$beta), c(.025, .500, 0.975))
  
  # It is quite close to the original CI, and the mean is correct
  
  ################################################################################
  #                                                                              #
  #  END                                                                         #
  #                                                                              #
  ################################################################################
  
  return(model)
  
}
  