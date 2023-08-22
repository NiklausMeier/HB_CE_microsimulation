################################################################################
#                                                                              #
# Study: Cost-Effectiveness Analysis of Gene Therapy for Haemophilia B         #
# Design: Cost-Effectiveness Model using bleeds as a continuous variable       #
# Outcome: Costs, QALYS, ICER                                                  #
# Task: Master file                                                            #
# Author: Niklaus Meier                                                        #
# R version: 4.2.1                                                             #
#                                                                              #
################################################################################

################################################################################
#                                                                              #
# STUDY SETUP                                                                  #
#                                                                              #
################################################################################

#===============================================================================
# Clear workspace
#===============================================================================

rm(list=ls())
save.image()
graphics.off()
gc()

#===============================================================================
# Load packages
#===============================================================================

packages <- c('ggplot2', 'RColorBrewer', 'pander', 'data.table', 'lookup',
              'truncnorm', 'ggpubr', 'plyr','ggrepel', 'captioner', 'rmarkdown')

invisible(lapply(packages, function(x)
  if( !require(x, character.only = TRUE)){
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  } else {
    library(x, character.only = TRUE)
  }))

rm(packages)

#===============================================================================
# Working directories
#===============================================================================

directories <- list()

directories$dir_dat_orig <-  './data_orig/' 
directories$dir_dat_deriv <-  './data_deriv/'

directories$dir_fun <- './functions/'
directories$dir_res <- './results/'

#===============================================================================
# Strings as factors
#===============================================================================

options(stringsAsFactors = FALSE)

#===============================================================================
# Turn off scientific notation
#===============================================================================

options(scipen = 999)

#===============================================================================
# Settings and global parameters
#===============================================================================

settings <- list() 

settings$seed <- 101010101

set.seed(settings$seed)

settings$timestamp <- format(Sys.time(), "%Y%m%d.%Hh%Mm%Ss")

settings$time <- list(yr2mon = 12.00,            # Months per year
                      yr2day = 365.25,           # Days per year
                      wk2day = 7,                # Days per week
                      yr2wk  = (365.25/7) ,      # Weeks per year
                      mon2wk = ((365.25/7)/12),  # Weeks per month
                      mon2day = (365.25/12)      # Days per month
)

#===============================================================================
# Ggplot theme
#===============================================================================

settings$ggplot_theme <- theme(panel.background = element_rect(fill = "white",
                                                               colour = "white",
                                                               linewidth = 0.5, 
                                                               linetype = "solid"),
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.grid.major.y = element_line(linewidth = 0.25, linetype = 'dashed', colour = "grey"),
                               text = element_text(size = 12),
                               plot.margin = margin(10, 10, 0, 10),
                               axis.line = element_line(linewidth = 0.25),
                               axis.text = element_text(size = 12),
                               #legend.title=element_text(size=12),
                               legend.title = element_blank(),
                               legend.text=element_text(size = 12),
                               legend.position = 'bottom',
                               legend.direction = "vertical",
                               legend.key = element_blank())

#===============================================================================
# Functions
#===============================================================================

## Load functions to prepare model
source(paste0(directories$dir_fun, "fun_model_setup.R"))
source(paste0(directories$dir_fun, "fun_parameter_baseline.R"))
source(paste0(directories$dir_fun, "fun_param_sample.R"))
source(paste0(directories$dir_fun, "fun_multinorminv.R"))
source(paste0(directories$dir_fun, "fun_scenarios.R"))

## Load functions to run model
source(paste0(directories$dir_fun, "fun_gen_pop.R"))
source(paste0(directories$dir_fun, "fun_ETRANACOGENE_abr.R"))
source(paste0(directories$dir_fun, "fun_PROPHYLAXIS_abr.R"))
source(paste0(directories$dir_fun, "fun_bleeds_death.R"))
source(paste0(directories$dir_fun, "fun_resources.R"))
source(paste0(directories$dir_fun, "fun_utility.R"))
source(paste0(directories$dir_fun, "fun_utility_scen_Grochtendreis.R"))
source(paste0(directories$dir_fun, "fun_costs.R"))

## Load functions to analyze results
source(paste0(directories$dir_fun, "fun_results.R"))
source(paste0(directories$dir_fun, "fun_diagnostics.R"))
source(paste0(directories$dir_fun, "fun_diagnostics_plots.R"))
source(paste0(directories$dir_fun, "fun_tornado.R"))

################################################################################
#                                                                              #
# RUN SCRIPTS                                                                  #
#                                                                              #
################################################################################

source('01_life_tables_weight.R', encoding = 'utf-8')

set.seed(settings$seed)
source('02A_cohort.R', encoding = 'utf-8')

rmarkdown::render('02B_analyze_cohort.Rmd',
                  output_file = paste(directories$dir_res,
                                      'cohort_',
                                      settings$timestamp,
                                      '.docx',
                                      sep = ''))

set.seed(settings$seed)
source('03A_deterministic_micro.R', encoding = 'utf-8')

rmarkdown::render('03B_analyze_deterministic_micro.Rmd',
                  output_file = paste(directories$dir_res,
                                      'deterministic_',
                                      settings$timestamp,
                                      '.docx',
                                      sep = ''))

set.seed(settings$seed)
source('04A_probabilistic_micro.R', encoding = 'utf-8')

rmarkdown::render('04B_analyze_probabilistic_micro.Rmd',
                  output_file = paste(directories$dir_res,
                                      'probabilistic_',
                                      settings$timestamp,
                                      '.docx',
                                      sep = ''))

set.seed(settings$seed)
source('05A_univariate_sensitivity_analysis.R', encoding = 'utf-8')

rmarkdown::render('05B_analyze_univariate_sensitivity_analysis.Rmd',
                  output_file = paste(directories$dir_res,
                                      'sens_analysis_',
                                      settings$timestamp,
                                      '.docx',
                                      sep = ''))

set.seed(settings$seed)
source('06A_scenarios.R', encoding = 'utf-8')

rmarkdown::render('06B_analyze_scenarios.Rmd',
                  output_file = paste(directories$dir_res,
                                      'scenarios_',
                                      settings$timestamp,
                                      '.docx',
                                      sep = ''))

