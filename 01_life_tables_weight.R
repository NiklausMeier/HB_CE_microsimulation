################################################################################
#                                                                              #
# Study: Cost-Effectiveness Analysis of Gene Therapy for Haemophilia B         #
# Design: Cost-Effectiveness Model using Bleeds as a continuous variable       #
# Outcome: Death Probabilities and Rates                                       #
# Task: Use German male life tables to calculate mortality based on age        #
# Author: Niklaus Meier                                                        #
# R version: 4.2.1                                                             #
#                                                                              #
################################################################################

start.time <- Sys.time()

################################################################################
#                                                                              #
# Setup                                                                        #
#                                                                              #
################################################################################

# Life Tables:
# Human Mortality Database. University of California, Berkeley (USA), and Max Planck Institute for Demographic Research (Germany). 
# Available at www.mortality.org (data downloaded on 26.01.2022). Germany, 1x1 Life Tables, Female and Male.

life_tables_male <- read.csv(paste0(directories$dir_dat_orig, '/Germany_life_tables_male.csv'), 
                             sep = ";", header = TRUE, fileEncoding = "UTF-8-BOM")

life_tables_female <- read.csv(paste0(directories$dir_dat_orig, '/Germany_life_tables_female.csv'), 
                               sep = ";", header = TRUE, fileEncoding = "UTF-8-BOM")
# Weight:
# Gesundheitsberichterstattung des Bundes (Germany, 2017). 
# Available at https://www.gbe-bund.de/gbe/!pkg_olap_tables.prc_set_orientation?p_uid=gast&p_aid=57250655&p_sprache=D&p_help=2&p_indnr=223&p_ansnr=26787077&p_version=4&D.000=3&D.002=2&D.003=1&D.100=1 (data downloaded on 26.01.2022).
weight <- read.csv(paste0(directories$dir_dat_orig, '/Germany_weight.csv'), 
                   sep = ";", header = TRUE, fileEncoding = "UTF-8-BOM")

################################################################################
#                                                                              #
# Life Tables Male                                                             #
#                                                                              #
################################################################################

# Only select data from the newest year, 2017.

life_tables_male <- life_tables_male[which(life_tables_male$Year=='2017'),]

# rename
names(life_tables_male)[names(life_tables_male) == "Age"] <- "age"
names(life_tables_male)[names(life_tables_male) == "qx"] <- "male_yearly_death_probability"

# Only keep these columns, get rid of rest

life_tables_male <- life_tables_male[c("age","male_yearly_death_probability")]

# Convert Age to numeric

life_tables_male$age <- as.numeric(life_tables_male$age)

# Calculate yearly death rate based on death probability

life_tables_male["male_yearly_death_rate"] <- -log(1-life_tables_male["male_yearly_death_probability"])

################################################################################
#                                                                              #
# Life Tables Female                                                           #
#                                                                              #
################################################################################

life_tables_female <- life_tables_female[which(life_tables_female$Year=='2017'),]

# rename
names(life_tables_female)[names(life_tables_female) == "Age"] <- "age"
names(life_tables_female)[names(life_tables_female) == "qx"] <- "female_yearly_death_probability"

# Only keep these columns, get rid of rest

life_tables_female <- life_tables_female[c("age","female_yearly_death_probability")]

# Convert Age to numeric

life_tables_female$age <- as.numeric(life_tables_female$age)

# Calculate yearly death rate based on death probability

life_tables_female["female_yearly_death_rate"] <- -log(1-life_tables_female["female_yearly_death_probability"])

################################################################################
#                                                                              #
# Life Tables and Weight Merged                                                #
#                                                                              #
################################################################################

life_tables <- merge(life_tables_male, life_tables_female, 
                     by = "age")

life_tables_weight <- merge(life_tables, weight, 
                            by = "age")

################################################################################
#                                                                              #
# Save Results                                                                 #
#                                                                              #
################################################################################

setDT(life_tables_weight)

save(life_tables_weight, file = paste0(directories$dir_dat_deriv, '/Germany_life_tables_weight.Rdata'))

################################################################################
#                                                                              #
# Cleanup                                                                      #
#                                                                              #
################################################################################

rm(life_tables_male)
rm(life_tables_female)
rm(weight)
rm(life_tables)
rm(life_tables_weight)

################################################################################
#                                                                              #
# REPORT RUNTIME                                                               #
#                                                                              #
################################################################################

end.time <- Sys.time()

print(paste('Run time =', round(as.numeric(end.time, units = "secs") - as.numeric(start.time, units = "secs"), 2)/60, 'minutes', sep = ' '))

gc()
