# Haemophilia B CE Microsimulation

## Description





## Installation
These files can be downladed and then run by a suitable version of R. The model was mainly developed in R 4.2.1 but should be compatible with future version, as long as the required packages do not change.  
The following packages (and their dependencies) are required: ggplot2, RColorBrewer, pander, data.table, lookup, truncnorm, ggpubr, plyr, ggrepel, captioner, rmarkdown  
These packages are installed and loaded in automatically when running the '00_master_file.R' script.  

## Usage

This section is intended as a brief overview for the operation of the R code used to perform the simulation of the HB cost-effectiveness model. It assumes prior familiarity with R. We will first describe the basic structure of the folder and their contents, followed by descriptions of the main scripts which execute the model, and conclude by describing the functions which are used by these scripts, and which contain most of the actual cost-effectiveness modelling in their code.

### Folders and files

The scripts and further files for this project are contained in a single folder, which has the following structure:  
•	**Main directory:** The main directory of the folder contains the scripts which execute the simulation, some additional files which R requires to operate, as well as folders which contain inputs or outputs of the model. We will briefly describe these folders top to bottom. 
	•	**data_deriv:** This folder contains outputs of the model as. Rdata files, which can again be loaded into R for further analysis.  
	•	**data_orig:** This folder contains inputs for the model from additional data sources, which can be loaded into the model.  
	•	**functions:** This folder contains the functions that the model requires.  
	•	**results:** This folder contains word files generated by the markdown scripts in the main directors.  

### R scripts and their usage

In the main folder of the project folder, there are multiple numbered scripts which can be run to execute the simulation. They are the following:  
•	**00_master_file.R:** This script should always be opened first, as it is required to setup the simulation as a whole. It deletes all previous data from the R workspace, loads in the required R packages, sets the directories  of the project, sets a theme for the generation of plots , loads in the required custom functions of the project, and determines various other settings. Under the final heading “RUN SCRIPTS” It can also be used to run all further scripts (01 – 06B) sequentially, though this can take a long time to perform, and is not required to run the individual further scripts.  
•	**01_life_tables_weight.R:** This script creates a table which contains yearly death probabilities and rates, as well as the average weight, for both men and women at each year of age from 0 to 110.  
•	**02A_cohort.R:** This script executes a simplified version of the model, in which a single average individual is simulated, rather than a heterogenous population. Therefore, it is essentially a cohort equivalent of the microsimulation. These results were not interpreted in the paper but seeing the outputs of a single individual is useful to check how the model functions and can be run much more rapidly than the microsimulation version.  
•	**02B_analyze_cohort.Rmd:** This script creates a word file analyzing the most important output of the cohort simulation.   
•	**03A_deterministic_micro.R:** This script executes a deterministic (or “base case”) version of the model, in which a heterogeneous population is simulated, but the other parameters of the model are fixed, rather than being probabilistically drawn.   
•	**03B_analyze_deterministic_micro.Rmd:** This script creates a word file analyzing the most important output of the deterministic simulation.  
•	**04A_probabilistic_micro:** This script executes a deterministic (or “base case”) version of the model, in which a heterogeneous population is simulated, and the other parameters of the model are probabilistically drawn.  
•	**04B_analyze_probabilistic_micro.Rmd:** This script creates a word file analyzing the most important output of the probabilistic simulation.  
•	**05A_univariate_sensitivity_analysis.R:** This script executes the deterministic model while altering individual parameters based on their probabilistic distributions.  
•	**05B_analyze_ univariate_sensitivity_analysis.Rmd:** This script creates a word file analyzing the most important output of the univariate sensitivity analysis.  
•	**06A_scenarios.R:** This script executes the deterministic model while altering individual parameters in a range and at intervals that are defined by the user to evaluate various scenarios.  
•	**06B_analyze_scenarios.Rmd:** This script creates a word file analyzing the most important output of the scenario analysis.  
These scripts, rather than duplicating the code for the model across the scripts, run a set of (mostly) identical functions, which contain the actual code that runs the model. This means that if the functions are altered, all of the scripts for the running of the model will be running this altered function, and the code does not need to be modified in multiple places.  

### Functions

The functions are loaded into R in the 00_master_file.R script. Most of these functions should be run sequentially to perform the model. Each of these functions performs a particular task, which is described below.  

* *Functions to prepare model:* * these functions must be run first as a preliminary step to create the necessary data structures, introduce data into the model, sample parameters, and define assumptions.  
•	**fun_model_setup.R:** This function requires the user to determine the cycle length, the number of cycles, the number of patients, the number of simulations (for probabilistic analysis), the treatments, and the willingness-to-pay threshold. Based on these inputs, the function creates a list to serve as the container for the relevant data, including the 3-dimensional arrays required to hold the patient data in the simulation. It also defines which measures (the second dimension of the array) are tracked by the model over time. If the user wishes to edit or add measures to the model, it should be done via this function, though changing the names of measures can cause errors if these names are not also adjusted in subsequent functions.  
•	**fun_parameter_baseline.R:** This function introduces the parameter values from assumptions and literature into the model, including base case values and further variables (standard deviations, standard errors, sample sizes, etc.) needed for the probabilistic simulation. The source of this data is also cited in the function. If the user wishes to edit or add further parameter values to the model, it should be done via this function.  
•	**fun_param_sample.R:** This function defines the statistical distributions for the probabilistic simulation and draws random values from these statistical distributions, based on the number of simulations which the user defined when setting up the model. If the user wishes to edit or add to the sampling from these statistical distributions, it should be done via this function.  
•	**fun_multinorminv.R:** This very simple function is used to randomly draw coefficients based on Cholesky decomposition of covariance matrix and is only used for the sampling of parameters.  
•	**fun_scenarios.R:** This function can be used to toggle certain modelling assumptions on or off, such as vial sharing or the ABR at which patients switch treatments.  

* *Functions to run model:* * these functions contain the actual model, including the population, the treatment strategies, the resulting bleeds and mortality outcomes, the resources consumption, health utilities, costs, and discounting. They must be run in the particular order in which they are listed to work properly.  
•	**fun_gen_pop.R:** This function generates a heterogeneous population which is then subsequently treated with the different strategies in subsequent functions. For each random parameter, a value is drawn from the statistical distribution to determine the value for each individual. This population size is automatically set to the size determined by the user when setting up the model. If the user simply wishes for the patients in the simulation to have baseline rather than randomly drawn parameter values, they can set the mode of the function to “homogeneous”. In practice, this is only used in the “02A_cohort.R” script.  
•	**Treatments:** Each of the treatment strategies has its own separate function, which determines the how and when treatments are administered as well as potentially switched. By having each treatment strategy in a separate function, it becomes easier to add or remove treatment strategies from the model at a later point in time. Based on the treatment strategy, the ABR is calculated in each cycle.  
**fun_ETRANACOGENE_abr.R:** to use this function, one of the treatment strategies set during model setup must be named “ETRANACOGENE”. In the first cycle, success or failure of the treatment with ED is randomly drawn and assigned to each patient. For those patients where treatment fails, they immediately switch treatment to prophylaxis for the rest of the model. For those patients where treatment succeeds, the ABR is calculated in all future cycles. In the first cycle where the ABR exceeds the treatment switch ABR threshold is exceeded, they switch treatment to FIX prophylaxis. ABR is then recalculated in subsequent cycles based on a multiplicative relative bleed reduction ED and FIX prophylaxis.  
**fun_PROPHYLAXIS_abr.R:** To use this function, one of the treatment strategies set during model setup must be named “PROPHYLAXIS”. Patients in the FIX prophylaxis strategy receive the same treatment in each cycle for the entirety of the model, such that their relative bleed reduction and thus ABR is constant across all cycles.  
•	**fun_bleeds_death.R:** This function calculates the number of bleeds in each cycle based on the ABR determined by the treatment strategy and the cycle length. These bleeds are either joint bleeds or other bleeds. The cumulative number of joint bleeds is calculated by adding the number of new joint bleeds to sum of the previous cycle. The Pettersson score is calculated based on the number of cumulative joint bleeds in each cycle. The patient-specific standardized mortality ratio is calculated by comparing the current ABR of each patient with the historic ABR and calculating a ratio. The patient-specific standardized mortality ratio is then multiplied by the background death rate (based on age) and converted into a probability, to receive the patient-specific probability of death in each cycle. The user can determine whether death should occur as an expected value or be randomly drawn in each cycle. If death is randomly drawn, then each living patient simply has a probability to die that cycle based on the above patient-specific probability of death. If death is an expected value, then the probability of being alive in each cycle is the product of multiplying the probability of being alive in the previous cycle by the probability of death in the previous cycle. Since patients can die at any point during a cycle, not just at the end of it, we perform a half-cycle correction by taking the average of the probability of being alive between each cycle and the one that follows it. The number of life years are calculated as the sum of the probability of being alive in each cycle.  
•	**fun_resources.R:** This function calculates the quantity of resources, in natural units, needed in each cycle. This encompasses the treatment with ED, prophylactic FIX, FIX for bleeds, hospitalizations for bleeds, and surgery based on the Pettersson score. Resource use is therefore determined by the history of each individual patient as well as the parameters of the simulation.  
•	**fun_utility.R:** This function calculates the quality of life of every patient in each cycle by summing up the baseline utility, additional utility for men, disutility for age and age-squared, disutility based on the Pettersson score, disutility for bleeds, disutility for surgery, and disutility if receiving factor infusion. This total utility is multiplied by the probability of being alive in each cycle to calculate the QALYs accumulated in that cycle.  
•	**fun_costs.R:** This function multiplies the resource use with the price per unit for each resource to calculate the cost due to that resource in each cycle.  
•	**fun_discounting.R:** This function multiplies undiscounted outcomes (life-years, QALYS, and costs) with the discount factor for that cycle to calculate the discounted outcome.  

* *Functions to analyze results:* * 
•	**fun_results.R:** this function summarizes and saves the results of the model into a separate data frame. It also calculates incremental QALYs, costs, NMB, and ICER. It can be run in two modes: “simulation” and “patient”.  
  **Simulation:** This mode saves the results of each simulation, aggregated from all patients in that simulation. This is especially important for the probabilistic analysis due to memory constraints when simulating a large number of patients.  
  **Patient:** This mode saves the results of each patient, across all simulations. For this to work properly, the function must be run after each individual simulation.  
•	**fun_diagnostics.R:** This function tracks the cumulative means of outcomes (across patients or simulations) and checks whether their values are reaching a stable convergence.  
•	**fun_diagnostics_plots.R:** This function generates plots for the cumulative means from the diagnostics to check for convergence visually.  
•	**fun_tornado.R:** This function generates a tornado plot based on a univariate sensitivity analysis from the “05A_univariate_sensitivity_analysis.R” script. For this, the user must input the results of the sensitivity analysis, the outcome which should be analyzed, and the baseline values of the results.  

An example of the R code utilizing these functions in order for the deterministic analysis would look like this:  

model <- fun_model_setup(cyc2day = settings$time$yr2day/4,  
                         ncycle = (92*4)+1,  
                         npatients = 10000,  
                         nsim = 1,  
                         mode = "deterministic",  
                         treatments = c(ETRANACOGENE", "PROPHYLAXIS"),  
                         threshold = 50000)  
model <- fun_parameter_baseline(model = model)  
model <- fun_param_sample(model = model)  
model <- fun_scenarios(model = model)  
model <- fun_gen_pop(model = model, mode = "heterogeneous")  
model <- fun_ETRANACOGENE_abr(model = model)  
model <- fun_PROPHYLAXIS_abr(model = model)  
model <- fun_bleeds_death(model = model, mode = "expected_value")  
model <- fun_resources(model = model)  
model <- fun_utility(model = model)  
model <- fun_costs(model = model)  
model <- fun_discounting(model = model)  
baseline_patient_results <- fun_results(model, mode = "patient")  
baseline_model_results  <- fun_results(model, mode = "simulation")  


## Support
If you have issues executing the code, or have other questions, please contact me at niklaus.meier@unibas.ch or 

## Authors and acknowledgment

Coder: Niklaus Meier, 
Code Review: Katya Galactionova
Further contributions to underlying theoretical work: Hendrik Fuchs, Cedric Hermans, Mark Pletscher, Matthias Schwenkglenks

## License
For open source projects, say how it is licensed.

## Project status
This model is currently complete, but may receive updates and further functionalities in the future.
