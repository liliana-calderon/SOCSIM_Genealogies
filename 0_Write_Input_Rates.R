#------------------------------------------------------------------------------------------------------
# SOCSIM - rsocsim rates - Write SOCSIM rate files from HFD and HMD data and presimulation .opop and .omar
# U:/SOCSIM/SOCSIM_Genealogies/0_Write_Input_Rates.R

## Write the input fertility and mortality rates for a SOCSIM micro-simulation for Sweden,
# using data from HFD (1891-2021) and HMD (1751-2021), with the HMDHFDplus package.
# and the Human Fertility Collection (HFC) for the period not covered in HFD
# Create a initial population file and empty marriage file for the simulations

# Created by Liliana Calderon on 08-06-2022
# Last modified by Liliana Calderon on 24-05-2023

#----------------------------------------------------------------------------------------------------
# Rate files format (Cf. Socsim oversimplified, p. 26):

# A rate block is a complete set of age speciﬁc rates governing a demographic event 
# for people of a particular sex, group and marital status. 
# In the rate files, the order matters and is always: 
# Event (birth, death, marriage), group (1, ..), sex (F or M), and
# marital status (married, single, divorced, widowed)
# For birth rates, this can be followed by number indicating parity.
# Each subsequent line contains a one month rate or probability
# and the age interval over which it holds (years and months of the upper age bound)
# The interval includes upper age bound in previous line and ends just before upper age bound in current line.
# The first two numbers are years and months of the upper age bound, which are added together.
# e.g. upper age bound of "0 1" refers to the first month of life, thus, [0,1)

#----------------------------------------------------------------------------------------------------
## Global Settings ----

# Clear work space
rm(list=ls(all=TRUE))

# Prevent scientific notation
options(scipen=999999)

# Load necessary packages
library(tidyverse)
library(readr)
# install.packages("HMDHFDplus")
library(HMDHFDplus)

# Load functions to write SOCSIM rate files from HFD/HMD and HFC
source("Functions_Input_Rates.R")

#----------------------------------------------------------------------------------------------------
## Write fertility rate files for SOCSIM using data from HFD ----

## NOTES:
# We keep the whole age range included in HFD [12-55], 
# but limit the open-ended age intervals 12- and 55+ to one year, i.e. [12-13) and [55-56)
# Additionally, SOCSIM requires a final age category with an upper bound
# corresponding to the maximum length of life. Here, we fixed it to [56-110]

# To convert the annual rates to monthly rates
# we assume that the fertility rates are constant over the year and divide them by 12 (months)
# According to SOCSIM oversimplified, these are rates rather than probabilities. 
# So multiplying a rate by the number of months in the age category gives
# the expected number of births that a woman who lives through the age category will experience.

# The rates are identical for married, single, divorced and widowed females. 
# We only specify the rates for single and married, 
# as other marital status will follow the rate default rules. 

# To write the input fertility rate files from HFD, 
# please type inside the quotes the name of the country as used in HFD) 
# and your HFD credentials (username and password) for the new website
# If needed, check the countries' names and availability with getHFDcountries()

write_socsim_rates_HFD(Country = "Type_here_Country_name", # SWE
                       HFD_username = "Type_here_HFD_username",
                       HFD_password = "Type_here_HFD_password")


#----------------------------------------------------------------------------------------------------
## Write fertility rate files for SOCSIM using data from Human Fertility Collection ----

# For the period before 1751-1890 that is covered by HMD but not by HFD, 
# We use data from the Human Fertility Collection
# Vital registers in Sweden prior 1891 came from Statistics Sweden (1969). 
# Historisk statistik för Sverige. Del 1. Bevolkning 1720-1967. Örebro : Statistiska centralbyrån. 
# we use Collection = "STAT" (SourceType = "Vital") 
# These data have Age Definition = ACY, Age in Completed Years and age range [15-49] 
# Age Interval is 1 for most ages [15-49], but for Age=14 (AgeInt=-99) and Age=55 (AgeInt=99)
# The age-specific rates hold over a period of 5 calendar years. 
# Hence, we use the same set of rates for each single year of the period.

write_socsim_rates_HFC(Country = "SWE", 
                       Collection = "STAT", Year_Max = 1891) # Vital statistics prior 1891 

#----------------------------------------------------------------------------------------------------
## Write mortality rate files for SOCSIM using data from HMD ----

## NOTES:

# We keep the whole age range included in HMD [0-110+]
# but limit the open-ended age interval 110+ to one year, i.e. [110-111)
# Originally, Socsim had an upper bound of 100 years, but the limit has been extended in rsocsim 

# As SOCSIM inputs, we use the probabilities of death (qx) from period life tables, 
# which are smoothed for ages 80 and above (cf. HMD explanatory notes on old-age mortality).
## To convert the annual probabilities to monthly probabilities,  
# we assume that the probability of dying is constant over the year
# and use the formula 1-(1-nqx)^(1/n) proposed by Kenneth Watcher (2014, p. 53). 
# For the open-ended age interval 110+ monthly probabilities are equal to 1/12

# The probabilities are identical for all marital status 
# (married, single, divorced and widowed) of each sex
# We only specify the rates by sex for single as other marital status will follow the rate default rules. 

# As explained by the MD explanatory note https://www.mortality.org/Data/ExplanatoryNotes 
# there are some difficulties in the calculation of mortality rates for ages 80+
# For period life tables, the central death rate m(x) is used to compute
# probabilities of death q(x).
# The values of m(x) below age 80 are by definition equal
# to the observed population death rate M(x) shown on each country page.
# At older ages, however, the number of deaths and the exposure-to-risk
# eventually become quite small,
# and thus observed death rates display considerable random variation.
# Therefore, we smooth the M(x) values for ages 80 and older
# and use these smoothed values to compute q(x) above a certain age
# (based on the number of observed deaths).
# For details, see the Methods Protocol (pp. 35-37).
# This procedure helps to avoid certain difficulties
# in period life table calculations at older ages that may be caused by:
# 1) extremely high death rates resulting from exposure being smaller than the number of deaths,
# 2) death rates of zero resulting from no deaths at an age where exposure is non-zero,
# 3) undefined death rates at all ages where exposure is zero.

# To write the input mortality rate files from HMD, 
# please type inside the quotes the name of the country (as used in HMD) 
# and your HMD credentials (username and password) for the new website
# If needed, check the countries' names and availability with getHMDcountries()

write_socsim_rates_HMD(Country = "Type_here_Country_name", # SWE
                       HMD_username = "Type_here_HMD_username",
                       HMD_password = "Type_here_HMD_password")

#------------------------------------------------------------------------------------------------------
# Create initial .opop and empty omar file for the simulation ----

# Set size of initial population
size_opop <-  50000

# Create data.frame with 14 columns and nrows = size_opop
presim.opop <- setNames(data.frame(matrix(data = 0, ncol = 14, nrow = size_opop)), 
                   c("pid","fem","group","nev","dob","mom","pop",
                     "nesibm","nesibp","lborn","marid","mstat","dod","fmult"))

# Add pid 1:sizeopop
presim.opop$pid <- 1:size_opop

# Add sex randomly
presim.opop$fem <- sample(0:1, nrow(presim.opop), replace = T)

# Add group 1 for all individuals
presim.opop$group <- 1

# Add random dates of birth (max age around 50)
presim.opop$dob <- sample(600:1200, nrow(presim.opop), replace = T)

## Check distribution of sex 
table(presim.opop$fem)

# Write initial population for pre-simulation (without fertility multiplier)
write.table(presim.opop, "presim.opop", row.names = F, col.names = F)


## Create an empty data frame for presim.omar
presim.omar <- data.frame()

# Write empty omar for pre-simulation
write.table(presim.omar, "presim.omar", row.names = F, col.names = F)