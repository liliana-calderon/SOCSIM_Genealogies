#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Genealogies - Run simulations 
# U:/SOCSIM/SOCSIM_Genealogies/1_Run_Simulations.R

## Run 10 SOCSIM demographic microsimulation using the 'rsocsim' package' and data for Sweden (1751-2021)
## and read the output into R. 

# Created by Liliana Calderon on 18-01-2022
# Last modified by Liliana Calderon on 26-05-2023

# NB: Some functions are based on external code and repositories specified under each section.
#------------------------------------------------------------------------------------------------------
# SOCSIM simulation: rate files and assumptions ----

# The SOCSIM simulation is implemented using the 'rsocsim' package,
# Available at https://github.com/MPIDR/rsocsim 

# It uses as input historical data for Sweden (1751-2022). 
# Input age-specific fertility rates come from the Human Fertility Collection (HFC) for 1751-1890,
# the Human Fertility Database (HFD)for the period 1891-2022
# and age-specific mortality rates (1751-2021) come from the Human Mortality Database.
# To run the simulation, original HFD and HMD rates are converted to monthly rates/probabilities
# and SOCSIM format using the 0_Write_Input_Rates.R script.

# Female fertility rates are identical for all marital status, but are specified for single and married women
# Other marital status (divorced, widowed, cohabiting) follow the SOCSIM rate default rules. 
# Mortality rates correspond to probabilities of death (qx) of HMD period life tables. 
# They are identical for all marital status of each sex, and are only specified for single women and men. 
# Other marital status will follow the rate default rules. 

# The first segment of the simulation runs for 100 years to produce a stable age structure, 
# based on 1751-HMD and 1751-HFC age-specific rates

#------------------------------------------------------------------------------------------------------
## Run SOCSIM simulations with 'rsocsim' package ----
# The instructions to install the rsocsim package and run the simulation 
# can be found on https://github.com/MPIDR/rsocsim/blob/main/readme.md
# To get the latest version from source, RTools and devtools must be installed.

# Clear work space
rm(list=ls(all=TRUE))

# library(devtools)

# Install rsocsim from Github with devtools:
# devtools::install_github("MPIDR/rsocsim")

# To get the updates, uninstall and install again the package
# remove.packages("rsocsim", lib="C:/Users/calderonbernal/AppData/Local/R/win-library/4.2")

# Load rsocsim package
library("rsocsim")

# Specify the working directory, where the supfile and ratefiles are.
# If the R session is running through the project file (SWE_example.Rproj), the following command can be used.
folder <- getwd()
# Otherwise, the working directory must be specified after "folder <- "

# Name of the supervisory file stored in the above folder:
supfile <- "socsim.sup"
# Sup file for rates retrieved from HFC/HFD and HMD (1751-2022), created with the 0_Write_Input_Rates.R,  
# using marry after childbirth and heterogeneous fertility directives

# Random number generator seed:
sims_seeds <- as.character(sample(1:99999, 10, replace = F))
# seed <- as.character(sample(1:99999,1)) #15829

## Run the simulations for the random seeds. 
for(seed in sims_seeds) {
  
  ### Run a single SOCSIM-simulation with a given folder and supervisory file,
  # using the "clustercall" process method, which allow to run several simulations. 
  rsocsim::socsim(folder, supfile, seed, process_method = "future")
  
}

# Save the seeds numbers to use them later to read the data
save(sims_seeds, file = "sims_seeds.Rda")
#----------------------------------------------------------------------------------------------------
## Read the output .opop file ----

# Load packages 
library(tidyverse)

# Load seeds numbers and convert them into numeric
load("sims_seeds.Rda")
sims_seeds <- as.numeric(sims_seeds)

# Iterate the function for the 10 seeds to read opop of the 10 simulations
sims_opop <- map(sims_seeds, ~ rsocsim::read_opop(folder = getwd(),
                                                  supfile = "socsim.sup",
                                                  seed = .,
                                                  suffix = "",
                                                  fn = NULL))

# Check the structure of sims_opop. A list of opop dfs
sims_opop %>% str()

# Save sims_opop list to use later
save(sims_opop, file = "sims_opop.RData")