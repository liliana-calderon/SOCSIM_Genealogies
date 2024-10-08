#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Genealogies - Run simulations 
# U:/SOCSIM/SOCSIM_Genealogies/1_Run_Simulations.R

## Run 10 SOCSIM demographic microsimulation using the 'rsocsim' package and data for Sweden (1751-2072)
## and read the output into R. 

# Created on 18-01-2022
# Last modified on 02-09-2024
#------------------------------------------------------------------------------------------------------
# SOCSIM simulation: rate files and assumptions ----

# The SOCSIM simulation is implemented using the 'rsocsim' package,
# Available at https://github.com/MPIDR/rsocsim 

# It uses as input historical data for Sweden (1751-2022). 
# Input age-specific fertility rates come from the Human Fertility Collection (HFC) for 1751-1890,
# the Human Fertility Database (HFD) for 1891-2022
# and age-specific mortality rates come from the Human Mortality Database for the whole period (1751-2022).
# To run the simulation, original HFC, HFD and HMD rates are converted to monthly rates/probabilities
# and SOCSIM format using the 0_Write_Input_Rates.R script.
# HFC fertility rates are provided by 5 calendar years and hence considered to be constant over each sub-period


# Female fertility rates are identical for all marital status, but are specified for single and married women
# Other marital status (divorced, widowed, cohabiting) follow the SOCSIM rate default rules. 
# Mortality rates correspond to probabilities of death (qx) of HMD period life tables. 
# They are identical for all marital status of each sex, and are only specified for single women and men. 
# Other marital status will follow the rate default rules. 

# The first segment of the simulation runs for 100 years to produce a stable age structure, 
# based on 1751-HFC and 1751-HMD age-specific rates
# The last segment of the simulation runs for another 100 years to avoid problems with the final population,
# based on 2022 HFD and HMD age-specific rates

#------------------------------------------------------------------------------------------------------
## Run SOCSIM simulations with 'rsocsim' package ----
# The instructions to install the rsocsim package and run the simulation 
# can be found on https://github.com/MPIDR/rsocsim/blob/main/readme.md
# To get the latest version from source, RTools and devtools must be installed.

# library(devtools)
# Install rsocsim from Github with devtools:
# devtools::install_github("MPIDR/rsocsim")

# To get the updates, uninstall and install again the package
# remove.packages("rsocsim")

# Load rsocsim package
library("rsocsim")

# Specify the working directory, where the supfile and ratefiles are.
# If the R session is running through the project file (SWE_example.Rproj), the following command can be used.
folder <- getwd()
# Otherwise, the working directory must be specified after "folder <- "

# Name of the supervisory file stored in the above folder:
supfile <- "Sweden_0_72.sup"
# Sup file for rates retrieved from HFC/HFD and HMD (1751-2072), created with the 0_Write_Input_Rates.R,  
# using marry after childbirth directive (heterogeneous fertility disabled)

# Random number generator seed:
sims_seeds <- as.character(sample(1:99999, 10, replace = F))

# NB: To get exactly the same simulation results, you should use the same randomly generated seeds
# sims_seeds <- c("51972", "16372", "47533", "10973", "52743", "23235", "75025", "46101", "29373", "76946")

# Save the seeds numbers to use them later to read the data
save(sims_seeds, file = "sims_seeds.Rda")

## Run the simulations for the random seeds. 
start <- Sys.time()
for(seed in sims_seeds) {
  
  ### Run a single SOCSIM-simulation with a given folder, seed and supervisory file
  rsocsim::socsim(folder, supfile, seed, process_method = "future")
  
}
end <- Sys.time()
print(end-start)
# Time difference of 4.39479 hours for 10 simulations (1651-2072), with initial population of 5000 and hetfert 0 
#----------------------------------------------------------------------------------------------------
## Read the output .opop and .omar files ----

# Load packages 
library(tidyverse)

# Load seeds numbers and convert them into numeric
load("sims_seeds.Rda")
sims_seeds <- as.numeric(sims_seeds)

# Iterate the function for the 10 seeds to read opop of the 10 simulations
sims_opop <- map(sims_seeds, ~ rsocsim::read_opop(folder = getwd(),
                                                  supfile = "Sweden_0_72.sup",
                                                  seed = .,
                                                  suffix = "",
                                                  fn = NULL))

# Check the structure of sims_opop. A list of opop dfs
sims_opop %>% str()

# Save sims_opop list to use later
save(sims_opop, file = "sims_opop.RData")


# Iterate the function for the 10 seeds to read omar of the 10 simulations
sims_omar <- map(sims_seeds, ~ rsocsim::read_omar(folder = getwd(),
                                                  supfile = "Sweden_0_72.sup",
                                                  seed = .,
                                                  suffix = "",
                                                  fn = NULL))

# Check the structure of sims_omar. A list of omar dfs
sims_omar %>% str()

# Save sims_opop list to use later
save(sims_omar, file = "sims_omar.RData")