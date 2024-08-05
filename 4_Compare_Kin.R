#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Sweden - Trace and compare genealogical subsets of kin of a given ego(s) 
# U:/SOCSIM/SOCSIM_Genealogies/4_Compare_Kin.R

## Trace relevant ascending and lateral kin up to the 4th generation
# of a given ego(s) from a SOCSIM microsimulation for Sweden (1751-2022)
# and compare demographic measures from the whole simulation and the genealogical subsets

# Created on 13-04-2022
# Last modified on 31-07-2024

## NB: To run this code, it is necessary to have already run the scripts 
# 1_Run_Simulations.R and 3_Compare_Ancestors.R
#------------------------------------------------------------------------------------------------------
## General settings and functions ----

# Prevent scientific notation (useful for the rate calculation)
options(scipen=999999)

## Load packages 
library(tidyverse)
library(ggh4x)  # To facet scales-
library(patchwork) # To combine ggplots
#library(rsocsim) # Functions to estimate rates
library(viridis)

## Load theme for the graphs and to convert SOCSIM time
source("Functions/Functions_Graphs.R")

## Load functions to estimate age-specific fertility and mortality rates 
# These are a slightly modified version of the functions in the rsocsim package 
# that allow to handle the intentional duplicates in the data (in fertility rates)
# and retrieve the last month from max(dod) instead of dob. 
# Important here as max(dob) in most subsets and simulations is not the last simulated month
# So, the function in rsocsim will assign incorrectly the last month of the simulation
# and hence calculate wrongly the years of birth and death
#source("Functions/Functions_Fertility_Rates_Mod.R")
source("Functions/Functions_Fertility_Rates_Mod_Denom.R") # Remove direct ancestors' offspring from denominator
source("Functions/Functions_Mortality_Rates_Mod.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with asmr calculated with estimate_mortality_rates_mod()
source("Functions/Functions_Life_Table.R")

## Load function to re-code the kin_type for the ancestors
source("Functions/Functions_Ancestors.R")

# Load modified function to retrieve kin
source("Functions/Functions_Kin_Mod.R")
#------------------------------------------------------------------------------------------------------
## Trace kin of people alive in 2023 ----

# Load saved list with opop from 10 simulations, generated in 1_Run_Simulations.R
load("sims_opop.RData")
# Load saved list with omar from 10 simulations, generated in 1_Run_Simulations.R
load("sims_omar.RData")

# Load same sample of egos for each opop used to trace the ancestors in 3_Compare_Ancestors.R
# These are a 10% sample of people alive at the end of the simulation, i.e. dod == 0, 
# who are older than 18 years old on 01-01-2023, i,e. dob 1914-2004 
load("Subsets/egos_samp_10.RData")

## Retrieve the relatives of each simulation sample of egos alive in 2023
kin_10_lc <- pmap_dfr(list(sims_opop, sims_omar, egos_samp_10),
                      ~ retrieve_kin_mod(opop = ..1, omar = ..2, pid = ..3,
                                         KidsOf = KidsOf),
                      .id = "Sim_id")
save(kin_10_lc, file = "Subsets/kin_10_lc.RData")
# Time difference of 3.9178 hours for 10 simulations with initial population of 5000, with hetfert0

# Select relevant variables from each opop and add Sim_id to allow merging with the kin_10 df 
sims_opop_id <- map_dfr(sims_opop, ~ select(.x, c(pid, fem, dob, dod, mom, lborn, marid, mstat)), 
                        .id = "Sim_id")

# Format the kin data,  put into long format and add columns from each opop
kin_10 <- kin_10_lc %>% 
  mutate(ego_id = unlist(egos_samp_10)) %>% 
  pivot_longer(-c(Sim_id, ego_id), names_to = "kin_type", values_to = "pid") %>% 
  unnest_longer(kin_type:pid) %>%
  filter(!is.na(pid) & pid != 0) %>% 
  # Sometimes, a pid is included as more than one kin type, 
  # e.g. gunclesaunts and ggparents of the same ego
  # Also,  "gggparents" and "ggunclesaunts", "ggggparents" and "gggunclesaunts" 
  # "gggggparents" and "ggggunclesaunts", "ggggggparents" and "gggggunclesaunts"
  # "gggggggparents" and "ggggggunclesaunts"
  # This come from strange kinship relations
  group_by(Sim_id, ego_id) %>% 
  distinct(pid, .keep_all = TRUE) %>% 
  ungroup() %>% 
  left_join(sims_opop_id, by = c("Sim_id", "pid"))

# Save the data frame with the relatives of 10 simulations samples
save(kin_10, file = "Subsets/kin_10.RData")

#------------------------------------------------------------------------------------------------------
## Merge direct ancestors with lateral kin of sample of egos alive in 2023 ----

# Load the data frame with the relatives of 10 simulations samples. 
# This has no duplicates from the same Sim_id and ego_id, 
# but still can have duplicates from different egos
load("Subsets/kin_10.RData")

# Load the data frame with the ancestors of the 10 simulations samples of egos alive in 2023
load("Subsets/ancestors_10.RData")

## Re-code the kin_type for the ancestors and remove duplicates (from the same simulation and ego)
ancestors_10 <- ancestors_10 %>% 
  recode_ancestors() %>%
  group_by(Sim_id, ego_id) %>% 
  distinct(pid, .keep_all = TRUE) %>% 
  ungroup() 

# Merge the two data frames and remove kin types included in both (i.e.g, parents, gparents, ggparents)
# NB: This data still contains duplicates within each simulation, as someone can be a relative of different egos
# Yet we will remove the duplicates only when creating the lists to estimate the rates
anc_kin_10 <- bind_rows(ancestors_10, kin_10) %>% 
  group_by(Sim_id, ego_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup()

# Save the data frame with the ancestors and relatives of 10 simulations samples
save(anc_kin_10, file = "Subsets/anc_kin_10.RData")

# Proportion by kin type, without duplicates by simulation
anc_kin_10 %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  pull(kin_type) %>% 
  table() %>% 
  prop.table()*100

#------------------------------------------------------------------------------------------------------
## Create lists with opop of selected kin keeping only unique pids ----

# Load the data frame with the ancestors and relatives of 10 simulations samples
load("Subsets/anc_kin_10.RData")

# All direct ancestors and their offspring:
# type_anc_off <- c("ego", "parents", "gparents", "ggparents", 
#                   "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
#                   "siblings", "unclesaunts", "gunclesaunts", "ggunclesaunts", 
#                   "gggunclesaunts", "ggggunclesaunts",  "gggggunclesaunts", "ggggggunclesaunts")

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_off <- anc_kin_10 %>%
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)


# Direct ancestors and siblings, 3.7%
type_anc_z <- c("ego", "parents", "gparents", "ggparents", 
                "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                "siblings") 

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_z <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_z) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)


# Direct ancestors, offspring above and aunts/uncles, 4.8%
type_anc_au <- c("ego", "parents", "gparents", "ggparents", 
                 "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                 "siblings", "unclesaunts")

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_au <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_au) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)


## Direct ancestors, offspring above and great-aunt/uncles, 6.3%
type_anc_gau <- c("ego", "parents", "gparents", "ggparents", 
                  "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                  "siblings", "unclesaunts", "gunclesaunts")

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_gau <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_gau) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)


# Direct ancestors, offspring above and great-great-aunt/uncles, 7.6%

type_anc_ggau <- c("ego", "parents", "gparents", "ggparents", 
                   "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                   "siblings", "unclesaunts", "gunclesaunts", "ggunclesaunts")

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_ggau <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_ggau) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)


# Direct ancestors, offspring above and great-great-great-aunt/uncles, 6.3 %
type_anc_gggau <- c("ego", "parents", "gparents", "ggparents", 
                    "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                    "siblings", "unclesaunts", "gunclesaunts", "ggunclesaunts", 
                    "gggunclesaunts")

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_gggau <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_gggau) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)


# Direct ancestors, offspring above and great-great-great-great-aunt/uncles, 4.5 %
type_anc_ggggau <- c("ego", "parents", "gparents", "ggparents", 
                     "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                     "siblings", "unclesaunts", "gunclesaunts", "ggunclesaunts", 
                     "gggunclesaunts", "ggggunclesaunts")

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_ggggau <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_ggggau) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)


# Direct ancestors, offspring above and great-great-great-great-great-aunt/uncles, 4.2 %
type_anc_gggggau <- c("ego", "parents", "gparents", "ggparents", 
                      "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                      "siblings", "unclesaunts", "gunclesaunts", "ggunclesaunts", 
                      "gggunclesaunts", "ggggunclesaunts", "gggggunclesaunts")

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_gggggau <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_gggggau) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)


# Direct ancestors, offspring above and great-great-great-great-great-aunt/uncles, 8.2 %
type_anc_ggggggau <- c("ego", "parents", "gparents", "ggparents", 
                       "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                       "siblings", "unclesaunts", "gunclesaunts", "ggunclesaunts", 
                       "gggunclesaunts", "ggggunclesaunts", "gggggunclesaunts", "ggggggunclesaunts")

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_ggggggau <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_ggggggau) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)

## Check maximum dob in the different subsets to check the last month is correctly retrieved. 

# All direct ancestors and their offspring:
anc_off %>% map(.f = ~ max(.x$dob)) %>% unique()
# Direct ancestors and siblings
anc_z %>% map(.f = ~ max(.x$dob)) %>% unique()
# Direct ancestors, offspring above and aunts/uncles
anc_au %>% map(.f = ~ max(.x$dob)) %>% unique()
## Direct ancestors, offspring above and great-aunt/uncles
anc_gau %>% map(.f = ~ max(.x$dob)) %>% unique()
# Direct ancestors, offspring above and great-great-aunt/uncles
anc_ggau %>% map(.f = ~ max(.x$dob)) %>% unique()
# Direct ancestors, offspring above and great-great-great-aunt/uncles
anc_gggau %>% map(.f = ~ max(.x$dob)) %>% unique()
# Direct ancestors, offspring above and great-great-great-great-aunt/uncles
anc_ggggau %>% map(.f = ~ max(.x$dob)) %>% unique()
# Direct ancestors, offspring above and great-great-great-great-great-aunt/uncles
anc_gggggau %>% map(.f = ~ max(.x$dob)) %>% unique()
# Direct ancestors, offspring above and great-great-great-great-great-aunt/uncles
anc_ggggggau %>% map(.f = ~ max(.x$dob)) %>% unique()

# The max(dob) in most subsets and simulations is not the last simulated month
# So, the function in rsocsim will assign incorrectly the last month of the simulation
# and hence calculate wrongly the years of birth and death
# We must use the modified rate functions (using dod) to ensure that the last month correctly retrieved
# With max(dod), in all subsets and simulations, the last month will be properly identified

#----------------------------------------------------------------------------------------------------
## Age-Specific Fertility and Mortality rates, 5x5  -----
# Retrieve and compare rates derived from the whole simulation with those from the genealogical subsets
# with direct ancestors and their offspring

##  Estimate ASFR and ASMR for the genealogical subsets with direct ancestors and their offspring

## All direct ancestors and their offspring

# Estimate age-specific fertility rates from the subset of all direct ancestors and their offspring
asfr_anc_off <- map_dfr(anc_off, ~ estimate_fertility_rates_mod(opop = .x,
                                                                final_sim_year = 2022, #[Jan-Dec]
                                                                year_min = 1750, # Closed [
                                                                year_max = 2020, # Open )
                                                                year_group = 5, 
                                                                age_min_fert = 10, # Closed [
                                                                age_max_fert = 55, # Open )
                                                                age_group = 5), # [,)
                        .id = "Sim_id") 
save(asfr_anc_off, file = "Measures/asfr_anc_off.RData")

# Estimate age-specific mortality rates from the subset of all direct ancestors and their offspring
asmr_anc_off <- map_dfr(anc_off, ~ estimate_mortality_rates_mod(opop = .x,
                                                                final_sim_year = 2022, #[Jan-Dec]
                                                                year_min = 1750, # Closed
                                                                year_max = 2020, # Open )
                                                                year_group = 5,
                                                                age_max_mort = 110, # Open )
                                                                age_group = 5), # [,)
                        .id = "Sim_id") 
save(asmr_anc_off, file = "Measures/asmr_anc_off.RData")


## Direct ancestors and siblings

# Estimate age-specific fertility rates from the subset of direct ancestors, offspring above and siblings
asfr_anc_z <- map_dfr(anc_z, ~ estimate_fertility_rates_mod(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed [
                                                            year_max = 2020, # Open )
                                                            year_group = 5, 
                                                            age_min_fert = 10, # Closed [
                                                            age_max_fert = 55, # Open )
                                                            age_group = 5), # [,)
                      .id = "Sim_id") 
save(asfr_anc_z, file = "Measures/asfr_anc_z.RData")

# Estimate age-specific mortality rates from the subset of direct ancestors, offspring above and siblings
asmr_anc_z <- map_dfr(anc_z, ~ estimate_mortality_rates_mod(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed
                                                            year_max = 2020, # Open )
                                                            year_group = 5,
                                                            age_max_mort = 110, # Open )
                                                            age_group = 5), # [,)
                      .id = "Sim_id") 
save(asmr_anc_z, file = "Measures/asmr_anc_z.RData")


## Direct ancestors, offspring above and aunts/uncles

# Estimate age-specific fertility rates from the subset of direct ancestors, offspring above and aunts/uncles
asfr_anc_au <- map_dfr(anc_au, ~ estimate_fertility_rates_mod(opop = .x,
                                                              final_sim_year = 2022, #[Jan-Dec]
                                                              year_min = 1750, # Closed [
                                                              year_max = 2020, # Open )
                                                              year_group = 5, 
                                                              age_min_fert = 10, # Closed [
                                                              age_max_fert = 55, # Open )
                                                              age_group = 5), # [,)
                       .id = "Sim_id") 
save(asfr_anc_au, file = "Measures/asfr_anc_au.RData")

# Estimate age-specific mortality rates from the subset of direct ancestors, offspring above and aunts/uncles
asmr_anc_au <- map_dfr(anc_au, ~ estimate_mortality_rates_mod(opop = .x,
                                                              final_sim_year = 2022, #[Jan-Dec]
                                                              year_min = 1750, # Closed
                                                              year_max = 2020, # Open )
                                                              year_group = 5,
                                                              age_max_mort = 110, # Open )
                                                              age_group = 5), # [,)
                       .id = "Sim_id") 
save(asmr_anc_au, file = "Measures/asmr_anc_au.RData")


## Direct ancestors, offspring above and great-aunt/uncles

# Estimate age-specific fertility rates from the subset of direct ancestors, offspring above and great-aunt/uncles
asfr_anc_gau <- map_dfr(anc_gau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                final_sim_year = 2022, #[Jan-Dec]
                                                                year_min = 1750, # Closed [
                                                                year_max = 2020, # Open )
                                                                year_group = 5, 
                                                                age_min_fert = 10, # Closed [
                                                                age_max_fert = 55, # Open )
                                                                age_group = 5), # [,)
                        .id = "Sim_id") 
save(asfr_anc_gau, file = "Measures/asfr_anc_gau.RData")

# Estimate age-specific mortality rates from the subset of direct ancestors, offspring above and great-aunt/uncles
asmr_anc_gau <- map_dfr(anc_gau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                final_sim_year = 2022, #[Jan-Dec]
                                                                year_min = 1750, # Closed
                                                                year_max = 2020, # Open )
                                                                year_group = 5,
                                                                age_max_mort = 110, # Open )
                                                                age_group = 5), # [,)
                        .id = "Sim_id") 
save(asmr_anc_gau, file = "Measures/asmr_anc_gau.RData")


## Direct ancestors, offspring above and great-great-aunt/uncles

# Estimate age-specific fertility rates from the subset of direct ancestors, offspring above and 2x-great-aunt/uncles
asfr_anc_ggau <- map_dfr(anc_ggau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                  final_sim_year = 2022, #[Jan-Dec]
                                                                  year_min = 1750, # Closed [
                                                                  year_max = 2020, # Open )
                                                                  year_group = 5, 
                                                                  age_min_fert = 10, # Closed [
                                                                  age_max_fert = 55, # Open )
                                                                  age_group = 5), # [,)
                         .id = "Sim_id") 
save(asfr_anc_ggau, file = "Measures/asfr_anc_ggau.RData")

# Estimate age-specific mortality rates from the subset of direct ancestors, offspring above and 2x-great-aunt/uncles
asmr_anc_ggau <- map_dfr(anc_ggau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                  final_sim_year = 2022, #[Jan-Dec]
                                                                  year_min = 1750, # Closed
                                                                  year_max = 2020, # Open )
                                                                  year_group = 5,
                                                                  age_max_mort = 110, # Open )
                                                                  age_group = 5), # [,)
                         .id = "Sim_id") 
save(asmr_anc_ggau, file = "Measures/asmr_anc_ggau.RData")


## Direct ancestors, offspring above and great-great-great-aunt/uncles

# Estimate age-specific fertility rates from the subset of direct ancestors, offspring above and 3x-great-great-aunt/uncles
asfr_anc_gggau <- map_dfr(anc_gggau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                    final_sim_year = 2022, #[Jan-Dec]
                                                                    year_min = 1750, # Closed [
                                                                    year_max = 2020, # Open )
                                                                    year_group = 5, 
                                                                    age_min_fert = 10, # Closed [
                                                                    age_max_fert = 55, # Open )
                                                                    age_group = 5), # [,)
                          .id = "Sim_id") 
save(asfr_anc_gggau, file = "Measures/asfr_anc_gggau.RData")

# Estimate age-specific mortality rates from the subset of direct ancestors, offspring above and 3x-great-great-aunt/uncles
asmr_anc_gggau <- map_dfr(anc_gggau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                    final_sim_year = 2022, #[Jan-Dec]
                                                                    year_min = 1750, # Closed
                                                                    year_max = 2020, # Open )
                                                                    year_group = 5,
                                                                    age_max_mort = 110, # Open )
                                                                    age_group = 5), # [,)
                          .id = "Sim_id") 
save(asmr_anc_gggau, file = "Measures/asmr_anc_gggau.RData")


## Direct ancestors, offspring above and great-great-great-great-aunt/uncles 

# Estimate age-specific fertility rates from the subset of direct ancestors, offspring above and 4x-great-aunt/uncles
asfr_anc_ggggau <- map_dfr(anc_ggggau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                      final_sim_year = 2022, #[Jan-Dec]
                                                                      year_min = 1750, # Closed [
                                                                      year_max = 2020, # Open )
                                                                      year_group = 5, 
                                                                      age_min_fert = 10, # Closed [
                                                                      age_max_fert = 55, # Open )
                                                                      age_group = 5), # [,)
                           .id = "Sim_id") 
save(asfr_anc_ggggau, file = "Measures/asfr_anc_ggggau.RData")

# Estimate age-specific mortality rates from the subset of direct ancestors, offspring above and 4x-great-aunt/uncles
asmr_anc_ggggau <- map_dfr(anc_ggggau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                      final_sim_year = 2022, #[Jan-Dec]
                                                                      year_min = 1750, # Closed
                                                                      year_max = 2020, # Open )
                                                                      year_group = 5,
                                                                      age_max_mort = 110, # Open )
                                                                      age_group = 5), # [,)
                           .id = "Sim_id") 
save(asmr_anc_ggggau, file = "Measures/asmr_anc_ggggau.RData")

# Direct ancestors, offspring above and great-great-great-great-great-aunt/uncles 

# Estimate age-specific fertility rates from the subset of direct ancestors, offspring above and 5x-great-aunt/uncles
asfr_anc_gggggau <- map_dfr(anc_gggggau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                        final_sim_year = 2022, #[Jan-Dec]
                                                                        year_min = 1750, # Closed [
                                                                        year_max = 2020, # Open )
                                                                        year_group = 5, 
                                                                        age_min_fert = 10, # Closed [
                                                                        age_max_fert = 55, # Open )
                                                                        age_group = 5), # [,)
                            .id = "Sim_id") 
save(asfr_anc_gggggau, file = "Measures/asfr_anc_gggggau.RData")

# Estimate age-specific mortality rates from the subset of direct ancestors, offspring above and 5x-great-aunt/uncles
asmr_anc_gggggau <- map_dfr(anc_gggggau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                        final_sim_year = 2022, #[Jan-Dec]
                                                                        year_min = 1750, # Closed
                                                                        year_max = 2020, # Open )
                                                                        year_group = 5,
                                                                        age_max_mort = 110, # Open )
                                                                        age_group = 5), # [,)
                                .id = "Sim_id") 
save(asmr_anc_gggggau, file = "Measures/asmr_anc_gggggau.RData")

## Direct ancestors, offspring above and great-great-great-great-great-aunt/uncles 

# Estimate age-specific fertility rates from the subset of direct ancestors, offspring above and 6x-great-aunt/uncles
asfr_anc_ggggggau <- map_dfr(anc_ggggggau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                          final_sim_year = 2022, #[Jan-Dec]
                                                                          year_min = 1750, # Closed [
                                                                          year_max = 2020, # Open )
                                                                          year_group = 5, 
                                                                          age_min_fert = 10, # Closed [
                                                                          age_max_fert = 55, # Open )
                                                                          age_group = 5), # [,)
                             .id = "Sim_id") 
save(asfr_anc_ggggggau, file = "Measures/asfr_anc_ggggggau.RData")

# Estimate age-specific mortality rates from the subset of direct ancestors, offspring above and 6x-great-aunt/uncles
asmr_anc_ggggggau <- map_dfr(anc_ggggggau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                          final_sim_year = 2022, #[Jan-Dec]
                                                                          year_min = 1750, # Closed
                                                                          year_max = 2020, # Open )
                                                                          year_group = 5,
                                                                          age_max_mort = 110, # Open )
                                                                          age_group = 5), # [,)
                             .id = "Sim_id") 
save(asmr_anc_ggggggau, file = "Measures/asmr_anc_ggggggau.RData")

#----------------------------------------------------------------------------------------------------
## Comparison of whole simulation with subsets of direct ancestors without and with their offspring ----

# ASFR ----

# Load ASFR 5x5 from the 10 simulations
load("Measures/asfr_10.RData")
# Load ASFR from the subset of "direct" family trees without duplicates
load("Measures/with_distinct/asfr_dir_wod.RData")
#load("Measures/asfr_dir_wod.RData")
# Load ASFR from the subset of all direct ancestors and their offspring
load("Measures/asfr_anc_off.RData")
# Load asfr from the subset of direct ancestors, offspring above and siblings
load("Measures/asfr_anc_z.RData")
# Load asfr from the subset of direct ancestors, offspring above and aunts/uncles
load("Measures/asfr_anc_au.RData")
# Load asfr from the subset of direct ancestors, offspring above and great-aunts/uncles
load("Measures/asfr_anc_gau.RData")
# Load asfr from the subset of direct ancestors, offspring above and 2x-great-aunts/uncles
load("Measures/asfr_anc_ggau.RData")
# Load asfr from the subset of direct ancestors, offspring above and 3x-great-aunts/uncles
load("Measures/asfr_anc_gggau.RData")
# Load asfr from the subset of direct ancestors, offspring above and 4x-great-aunts/uncles
load("Measures/asfr_anc_ggggau.RData")
# Load asfr from the subset of direct ancestors, offspring above and 5x-great-aunts/uncles
load("Measures/asfr_anc_gggggau.RData")
# Load asfr from the subset of direct ancestors, offspring above and 6x-great-aunts/uncles
load("Measures/asfr_anc_ggggggau.RData")

## Calculate the mean of the different simulations and add relevant columns

# Whole SOCSIM simulations
asfr_whole2 <- asfr_10 %>%
  group_by(year, age) %>% 
  summarise(ASFR = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Dataset = "Whole Simulation", 
         Rate = "ASFR") 

# Only direct ancestors without duplicates
asfr_dir_wod2 <- asfr_dir_wod %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct Ancestors",
         Rate = "ASFR") 

# All direct ancestors and their offspring
asfr_anc_off2 <- asfr_anc_off %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct Ancestors and their Offspring",
         Rate = "ASFR") 

# Direct ancestors and siblings
asfr_anc_z2 <- asfr_anc_z %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "1. + Siblings",
         Rate = "ASFR") 

# direct ancestors, offspring above and aunts/uncles
asfr_anc_au2 <- asfr_anc_au %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "2. + Aunts/Uncles",
         Rate = "ASFR") 

# Direct ancestors, offspring above and great-aunts/uncles
asfr_anc_gau2 <- asfr_anc_gau %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "3. + Great-Aunts/Uncles",
         Rate = "ASFR") 

# Direct ancestors, offspring above and 2x-great-aunts/uncles
asfr_anc_ggau2 <- asfr_anc_ggau %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "4. + 2x-Great-Aunts/Uncles",
         Rate = "ASFR") 

# Direct ancestors, offspring above and 3x-great-aunts/uncles
asfr_anc_gggau2 <- asfr_anc_gggau %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "5. + 3x-Great-Aunts/Uncles",
         Rate = "ASFR") 

# Direct ancestors, offspring above and 4x-great-aunts/uncles
asfr_anc_ggggau2 <- asfr_anc_ggggau %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "6. + 4x-Great-Aunts/Uncles",
         Rate = "ASFR") 

# Direct ancestors, offspring above and 5x-great-aunts/uncles
asfr_anc_gggggau2 <- asfr_anc_gggggau %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "7. + 5x-Great-Aunts/Uncles",
         Rate = "ASFR") 

# Direct ancestors, offspring above and 6x-great-aunts/uncles
asfr_anc_ggggggau2 <- asfr_anc_ggggggau %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "8. + 6x-Great-Aunts/Uncles",
         Rate = "ASFR") 


## Plot ASFR from whole SOCSIM simulation and subsets of direct and extended family trees without duplicates

# Same years to plot than above (in intervals). 
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)")

bind_rows(asfr_whole2, asfr_dir_wod2, asfr_anc_z2, asfr_anc_au2, asfr_anc_gau2, asfr_anc_ggau2, 
          asfr_anc_gggau2, asfr_anc_ggggau2, asfr_anc_gggggau2, asfr_anc_ggggggau2) %>%
  filter(year %in% yrs_plot & !is.nan(ASFR)) %>%
  mutate(Dataset = factor(Dataset, levels =  c("Direct Ancestors", "1. + Siblings", "2. + Aunts/Uncles", "3. + Great-Aunts/Uncles", 
                                               "4. + 2x-Great-Aunts/Uncles", "5. + 3x-Great-Aunts/Uncles", "6. + 4x-Great-Aunts/Uncles", 
                                               "7. + 5x-Great-Aunts/Uncles", "8. + 6x-Great-Aunts/Uncles", "Whole Simulation"))) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset), colour = Dataset))+
  facet_grid(year ~ . , scales = "free") +
  geom_line(linewidth = 1.2)+
  geom_point(size = 9)+
  scale_color_viridis(option = "H", discrete = T, direction = -1)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp2_ASFR.jpeg", width=17, height=14, dpi=300)

# ASMR ----

# Load ASMR 5x5 from the 10 simulations
load("Measures/asmr_10.RData")
# Load ASMR from the subset of "direct" family trees without duplicates
load("Measures/with_distinct/asmr_dir_wod.RData")
# Load asmr from the subset of all direct ancestors and their offspring
load("Measures/asmr_anc_off.RData")
# Load asmr from the subset of direct ancestors, offspring above and siblings
load("Measures/asmr_anc_z.RData")
# Load asmr from the subset of  direct ancestors, offspring above and aunts/uncles
load("Measures/asmr_anc_au.RData")
# Load asmr from the subset of direct ancestors, offspring above and great-aunts/uncles
load("Measures/asmr_anc_gau.RData")
# Load asmr from the subset of direct ancestors, offspring above and 2x-great-aunts/uncles
load("Measures/asmr_anc_ggau.RData")
# Load asmr from the subset of direct ancestors, offspring above and 3x-great-aunts/uncles
load("Measures/asmr_anc_gggau.RData")
# Load asmr from the subset of direct ancestors, offspring above and 4x-great-aunts/uncles
load("Measures/asmr_anc_ggggau.RData")
# Load asmr from the subset of direct ancestors, offspring above and 5x-great-aunts/uncles
load("Measures/asmr_anc_gggggau.RData")
# Load asmr from the subset of direct ancestors, offspring above and 6x-great-aunts/uncles
load("Measures/asmr_anc_ggggggau.RData")

# Whole SOCSIM simulations
asmr_whole2 <- asmr_10 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Whole Simulation",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# "Direct" family trees without duplicates
asmr_dir_wod2 <- asmr_dir_wod %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Direct Ancestors",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All direct ancestors and their offspring
asmr_anc_off2 <- asmr_anc_off %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "Direct Ancestors and their Offspring",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Direct ancestors and siblings
asmr_anc_z2 <- asmr_anc_z %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "1. + Siblings",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# direct ancestors, offspring above and aunts/uncles
asmr_anc_au2 <- asmr_anc_au %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "2. + Aunts/Uncles",
         Rate = "ASMR") %>%   
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# direct ancestors, offspring above and great-aunts/uncles
asmr_anc_gau2 <- asmr_anc_gau %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "3. + Great-Aunts/Uncles",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# direct ancestors, offspring above and 2x-great-aunts/uncles
asmr_anc_ggau2 <- asmr_anc_ggau %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "4. + 2x-Great-Aunts/Uncles",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# direct ancestors, offspring above and 3x-great-aunts/uncles
asmr_anc_gggau2 <- asmr_anc_gggau %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "5. + 3x-Great-Aunts/Uncles",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# direct ancestors, offspring above and 4x-great-aunts/uncles
asmr_anc_ggggau2 <- asmr_anc_ggggau %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "6. + 4x-Great-Aunts/Uncles",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# direct ancestors, offspring above and 5x-great-aunts/uncles
asmr_anc_gggggau2 <- asmr_anc_gggggau %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "7. + 5x-Great-Aunts/Uncles",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# direct ancestors, offspring above and 6x-great-aunts/uncles
asmr_anc_ggggggau2 <- asmr_anc_ggggggau %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "8. + 6x-Great-Aunts/Uncles",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

## Plotting ASMR from whole SOCSIM simulation and subsets of "direct" and "extended" family trees without duplicates

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_dir_wod$age)

# Same years to plot than above (in intervals). 
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asmr_whole2, asmr_dir_wod2, asmr_anc_z2, asmr_anc_au2, asmr_anc_gau2, asmr_anc_ggau2, 
          asmr_anc_gggau2, asmr_anc_ggggau2, asmr_anc_gggggau2, asmr_anc_ggggggau2)  %>%
  rename(Year = year) %>%
  filter(Year %in% yrs_plot) %>%
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>%
  mutate(Dataset = factor(Dataset, levels =  c("Direct Ancestors", "1. + Siblings", "2. + Aunts/Uncles", "3. + Great-Aunts/Uncles", 
                                               "4. + 2x-Great-Aunts/Uncles", "5. + 3x-Great-Aunts/Uncles", "6. + 4x-Great-Aunts/Uncles", 
                                               "7. + 5x-Great-Aunts/Uncles", "8. + 6x-Great-Aunts/Uncles", "Whole Simulation"))) %>%
  ggplot(aes(x = age, y = mx, group = interaction(Year, Dataset), colour = Dataset))+
  facet_grid(Year ~ Sex) +
  geom_line(linewidth = 1.2, show.legend = T)+
  geom_point(size = 9)+
  scale_color_viridis(option = "H", discrete = T, direction = -1)+
  theme_graphs()+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp2_ASMR.jpeg", width=17, height=15, dpi=300)

#----------------------------------------------------------------------------------------------------
## Final plot combining ASFR and ASMR ----

# Choose one year and age groups to plot
yrs_plot1 <- c("[1900,1905)") 
age_plot <- c("[0,1)", "[1,5)", "[10,15)", "[20,25)", "[30,35)", "[40,45)", "[50,55)",  "[60,65)", 
              "[70,75)", "[80,85)", "[90,95)", "[100,105)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_whole2$age)

# Define the same y breaks for all plots
y_breaks_asfr <- c(0.0, 0.05, 0.1, 0.15, 0.2)
y_breaks_asmr <- c(0.0, 0.0001, 0.001, 0.01, 0.1, 0.3)

## Plotting ASFR and ASMR (for females) from whole SOCSIM simulation and subsets of direct ancestors and their offspring
By_Age_Exp2 <- 
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR),
          asfr_dir_wod2 %>% rename(Estimate = ASFR),
          asfr_anc_off2 %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx),
            asmr_dir_wod2 %>% rename(Estimate = mx),
            asmr_anc_off2 %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot1) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels),
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates"), 
         Dataset = factor(Dataset, levels =  c("Direct Ancestors", "Direct Ancestors and their Offspring", "Whole Simulation"))) %>%
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(linewidth = 1.3, show.legend = T)+
  geom_point(data = . %>% filter(age %in% age_plot), 
             aes(shape = Dataset), size = 11) +
  scale_color_manual(values = c("#2779B7"))+ 
  scale_shape_manual(values = c(19, 18, 46)) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() + theme(legend.title = element_text(size = 18)) +
  labs(x = "Age") +
  guides(shape = guide_legend(order = 1), col = guide_legend(order = 2)) +
  theme(legend.justification = "left", 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
By_Age_Exp2

## Plot ASFR and ASMR (for females), with three years for appendix

# Choose three years to plot
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2 %>% rename(Estimate = ASFR),
          asfr_dir_wod2 %>% rename(Estimate = ASFR),
          asfr_anc_off2 %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx),
            asmr_dir_wod2 %>% rename(Estimate = mx),
            asmr_anc_off2 %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels),
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates"), 
         Dataset = factor(Dataset, levels =  c("Direct Ancestors", "Direct Ancestors and their Offspring", "Whole Simulation"))) %>%
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(Year ~ Rate, nrow = 3, ncol = 2, scales = "free") + 
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(data = . %>% filter(age %in% age_plot), 
             aes(shape = Dataset), size = 11) +
  scale_color_manual(values = c("#79B727", "#2779B7", "#B72779"))+ 
  scale_shape_manual(values = c(19, 18, 46)) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10"),
                               ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10"),
                               ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() + theme(legend.title = element_text(size = 18)) +
  labs(x = "Age") +
  guides(shape = guide_legend(order = 1), col = guide_legend(order = 2)) +
  theme(legend.justification = "left", 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
ggsave(file="Graphs/Socsim_Exp2_ASFR_ASMR_grouped.jpeg", width=18, height=25, dpi=300)

#----------------------------------------------------------------------------------------------------
# Figure for EPC presentation

# Choose three years to plot
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 
age_plot <- c("[0,1)", "[1,5)", "[10,15)", "[20,25)", "[30,35)", "[40,45)", "[50,55)",  "[60,65)", 
              "[70,75)", "[80,85)", "[90,95)", "[100,105)") 

bind_rows(asmr_whole2, asmr_dir_wod2, asmr_anc_off2) %>% 
  rename(Year = year) %>% 
  filter(Year %in% yrs_plot) %>%
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>% 
  ggplot(aes(x = age, y = mx, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap( ~ Year, nrow = 1, ncol = 3, scales = "free") +
  geom_line(linewidth = 1.3, show.legend = TRUE)+
  geom_point(data = . %>% filter(age %in% age_plot), 
             aes(shape = Dataset), size = 11) +
  scale_color_manual(values = c("#79B727", "#2779B7", "#B72779"))+ 
  scale_shape_manual(values = c(19, 18, 46)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs() +
  labs(x = "Age") +
  theme(legend.justification = "left", 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  guides(shape = guide_legend(order = 1))

ggsave(file="Graphs/Socsim_Exp2_ASMR_years.jpeg", width=24, height=9, dpi=300)
#----------------------------------------------------------------------------------------------------
## Summary measures: TFR and e0 ----
# Here, we use the rates by 1 year age group and 1 calendar year

# Estimate Total Fertility Rate from asfr 1x1 ----

# Estimate age-specific fertility rates 1x1 from the subset of all direct ancestors and their offspring
asfr_anc_off_1 <- map_dfr(anc_off, ~ estimate_fertility_rates_mod(opop = .x,
                                                                  final_sim_year = 2022, #[Jan-Dec]
                                                                  year_min = 1750, # Closed [
                                                                  year_max = 2023, # Open )
                                                                  year_group = 1, 
                                                                  age_min_fert = 10, # Closed [
                                                                  age_max_fert = 55, # Open )
                                                                  age_group = 1), # [,)
                          .id = "Sim_id") 
save(asfr_anc_off_1, file = "Measures/asfr_anc_off_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset of direct ancestors, offspring above and siblings
asfr_anc_z_1 <- map_dfr(anc_z, ~ estimate_fertility_rates_mod(opop = .x,
                                                              final_sim_year = 2022, #[Jan-Dec]
                                                              year_min = 1750, # Closed [
                                                              year_max = 2023, # Open )
                                                              year_group = 1, 
                                                              age_min_fert = 10, # Closed [
                                                              age_max_fert = 55, # Open )
                                                              age_group = 1), # [,)
                        .id = "Sim_id") 
save(asfr_anc_z_1, file = "Measures/asfr_anc_z_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset of direct ancestors, offspring above and aunts/uncles
asfr_anc_au_1 <- map_dfr(anc_au, ~ estimate_fertility_rates_mod(opop = .x,
                                                                final_sim_year = 2022, #[Jan-Dec]
                                                                year_min = 1750, # Closed [
                                                                year_max = 2023, # Open )
                                                                year_group = 1, 
                                                                age_min_fert = 10, # Closed [
                                                                age_max_fert = 55, # Open )
                                                                age_group = 1), # [,)
                         .id = "Sim_id") 
save(asfr_anc_au_1, file = "Measures/asfr_anc_au_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset of direct ancestors, offspring above and great-aunts/uncles
asfr_anc_gau_1 <- map_dfr(anc_gau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                  final_sim_year = 2022, #[Jan-Dec]
                                                                  year_min = 1750, # Closed [
                                                                  year_max = 2023, # Open )
                                                                  year_group = 1, 
                                                                  age_min_fert = 10, # Closed [
                                                                  age_max_fert = 55, # Open )
                                                                  age_group = 1), # [,)
                          .id = "Sim_id") 
save(asfr_anc_gau_1, file = "Measures/asfr_anc_gau_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset of direct ancestors, offspring above and 2x-great-aunts/uncles
asfr_anc_ggau_1 <- map_dfr(anc_ggau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                    final_sim_year = 2022, #[Jan-Dec]
                                                                    year_min = 1750, # Closed [
                                                                    year_max = 2023, # Open )
                                                                    year_group = 1, 
                                                                    age_min_fert = 10, # Closed [
                                                                    age_max_fert = 55, # Open )
                                                                    age_group = 1), # [,)
                          .id = "Sim_id") 
save(asfr_anc_ggau_1, file = "Measures/asfr_anc_ggau_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset of direct ancestors, offspring above and 3x-great-aunts/uncles
asfr_anc_gggau_1 <- map_dfr(anc_gggau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                      final_sim_year = 2022, #[Jan-Dec]
                                                                      year_min = 1750, # Closed [
                                                                      year_max = 2023, # Open )
                                                                      year_group = 1, 
                                                                      age_min_fert = 10, # Closed [
                                                                      age_max_fert = 55, # Open )
                                                                      age_group = 1), # [,)
                            .id = "Sim_id") 
save(asfr_anc_gggau_1, file = "Measures/asfr_anc_gggau_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset of direct ancestors, offspring above and 4x-great-aunts/uncles
asfr_anc_ggggau_1 <- map_dfr(anc_ggggau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                        final_sim_year = 2022, #[Jan-Dec]
                                                                        year_min = 1750, # Closed [
                                                                        year_max = 2023, # Open )
                                                                        year_group = 1, 
                                                                        age_min_fert = 10, # Closed [
                                                                        age_max_fert = 55, # Open )
                                                                        age_group = 1), # [,)
                             .id = "Sim_id") 
save(asfr_anc_ggggau_1, file = "Measures/asfr_anc_ggggau_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset of direct ancestors, offspring above and 5x-great-aunts/uncles
asfr_anc_gggggau_1 <- map_dfr(anc_gggggau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                          final_sim_year = 2022, #[Jan-Dec]
                                                                          year_min = 1750, # Closed [
                                                                          year_max = 2023, # Open )
                                                                          year_group = 1, 
                                                                          age_min_fert = 10, # Closed [
                                                                          age_max_fert = 55, # Open )
                                                                          age_group = 1), # [,)
                              .id = "Sim_id") 
save(asfr_anc_gggggau_1, file = "Measures/asfr_anc_gggggau_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset of direct ancestors, offspring above and 6x-great-aunts/uncles
asfr_anc_ggggggau_1 <- map_dfr(anc_ggggggau, ~ estimate_fertility_rates_mod(opop = .x,
                                                                            final_sim_year = 2022, #[Jan-Dec]
                                                                            year_min = 1750, # Closed [
                                                                            year_max = 2023, # Open )
                                                                            year_group = 1, 
                                                                            age_min_fert = 10, # Closed [
                                                                            age_max_fert = 55, # Open )
                                                                            age_group = 1), # [,)
                               .id = "Sim_id") 
save(asfr_anc_ggggggau_1, file = "Measures/asfr_anc_ggggggau_1.RData")

# Load ASFR 1x1 and calculate TFR for plotting ----

# Load asfr 1x1 from the 10 simulations
load("Measures/asfr_10_1.RData")
# Load asfr 1x1 from the subset of direct ancestors
load("Measures/with_distinct/asfr_dir_wod_1.RData")
# Load asfr 1x1 from the subset of direct ancestors and their offspring
load("Measures/asfr_anc_off_1.RData")
# Load asfr 1x1 from the subset of direct ancestors and siblings
load("Measures/asfr_anc_z_1.RData")
# Load asfr 1x1 from the subset of direct ancestors, offspring above and aunts/uncles
load("Measures/asfr_anc_au_1.RData")
# Load asfr 1x1 from the subset of direct ancestors, offspring above and great-aunts/uncles
load("Measures/asfr_anc_gau_1.RData")
# Load asfr 1x1 from the subset of direct ancestors, offspring above and 2x-great-aunts/uncles
load("Measures/asfr_anc_ggau_1.RData")
# Load asfr 1x1 from the subset of direct ancestors, offspring above and 3x-great-aunts/uncles
load("Measures/asfr_anc_gggau_1.RData")
# Load asfr 1x1 from the subset of direct ancestors, offspring above and 4x-great-aunts/uncles
load("Measures/asfr_anc_ggggau_1.RData")
# Load asfr 1x1 from the subset of direct ancestors, offspring above and 5x-great-aunts/uncles
load("Measures/asfr_anc_gggggau_1.RData")
# Load asfr 1x1 from the subset of direct ancestors, offspring above and 6x-great-aunts/uncles
load("Measures/asfr_anc_ggggggau_1.RData")

# Age breaks of fertility rates. Extract all the unique numbers from the intervals 
age_breaks_fert_1 <- unique(as.numeric(str_extract_all(asfr_10_1$age, "\\d+", simplify = T)))

# Retrieve age_group size
age_group_fert_1 <- unique(diff(age_breaks_fert_1))

# Whole SOCSIM simulations
TFR_whole <- asfr_10_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "Whole Simulation",
         Rate = "TFR", 
         sex = "female")

# Direct ancestors without duplicates
TFR_dir_wod <- asfr_dir_wod_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "Direct Ancestors",
         Rate = "TFR", 
         sex = "female")

# All direct ancestors and their offspring
TFR_anc_off <- asfr_anc_off_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "Direct Ancestors and their Offspring",
         Rate = "TFR",           
         sex = "female")

# Direct ancestors and siblings
TFR_anc_z <- asfr_anc_z_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "1. + Siblings",
         Rate = "TFR",           
         sex = "female") 

# direct ancestors, offspring above and aunts/uncles
TFR_anc_au <- asfr_anc_au_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "2. + Aunts/Uncles",
         Rate = "TFR",           
         sex = "female") 

# direct ancestors, offspring above and great-aunts/uncles
TFR_anc_gau <- asfr_anc_gau_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "3. + Great-Aunts/Uncles",
         Rate = "TFR",           
         sex = "female") 

# direct ancestors, offspring above and 2x-great-aunts/uncles
TFR_anc_ggau <- asfr_anc_ggau_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "4. + 2x-Great-Aunts/Uncles",
         Rate = "TFR",           
         sex = "female") 

# direct ancestors, offspring above and 3x-great-aunts/uncles
TFR_anc_gggau <- asfr_anc_gggau_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "5. + 3x-Great-Aunts/Uncles",
         Rate = "TFR",           
         sex = "female") 

# direct ancestors, offspring above and 4x-great-aunts/uncles
TFR_anc_ggggau <- asfr_anc_ggggau_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "6. + 4x-Great-Aunts/Uncles",
         Rate = "TFR",           
         sex = "female") 

# direct ancestors, offspring above and 5x-great-aunts/uncles
TFR_anc_gggggau <- asfr_anc_gggggau_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "7. + 5x-Great-Aunts/Uncles",
         Rate = "TFR",           
         sex = "female") 

# direct ancestors, offspring above and 6x-great-aunts/uncles
TFR_anc_ggggggau <- asfr_anc_ggggggau_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "8. + 6x-Great-Aunts/Uncles",
         Rate = "TFR",           
         sex = "female") 

## Plot TFR from whole SOCSIM simulation and subsets of "direct" and "extended"  family trees without duplicates
## Here we need to compute the mean for each dataset

bind_rows(TFR_whole, TFR_dir_wod  %>% filter(TFR != 0 & !is.nan(TFR)), 
          TFR_anc_z, TFR_anc_au, TFR_anc_gau, TFR_anc_ggau, 
          TFR_anc_gggau, TFR_anc_ggggau, TFR_anc_gggggau, TFR_anc_ggggggau) %>%
  mutate(Dataset = factor(Dataset, levels =  c("Direct Ancestors", "1. + Siblings", "2. + Aunts/Uncles", "3. + Great-Aunts/Uncles", 
                                               "4. + 2x-Great-Aunts/Uncles", "5. + 3x-Great-Aunts/Uncles", "6. + 4x-Great-Aunts/Uncles", 
                                               "7. + 5x-Great-Aunts/Uncles", "8. + 6x-Great-Aunts/Uncles", "Whole Simulation"))) %>%
  group_by(Year, Dataset) %>% 
  summarise(TFR = mean(TFR, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year, y = TFR)) +
  geom_line(aes(colour = Dataset), linewidth = 1.3)+
  scale_color_viridis(option = "H", discrete = T, direction = -1)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp2_TFR.jpeg", width=17, height=9, dpi=300)

# Summary measure of error in TFR ----

# Difference in means
DiM_TFR_Exp2 <- bind_rows(TFR_whole, TFR_dir_wod, TFR_anc_z, TFR_anc_au, TFR_anc_gau, TFR_anc_ggau, 
                          TFR_anc_gggau, TFR_anc_ggggau, TFR_anc_gggggau, TFR_anc_ggggggau) %>%
  filter(Year > 1750) %>% 
  # There can be rates of 0 and NaN values as the genealogist (most recent generation) are at least 18 years old. 
  # This happens after 2004
  filter(TFR != 0 & !is.nan(TFR)) %>% 
  group_by(Year, Dataset) %>% 
  summarise(TFR = mean(TFR, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = Year, names_from = "Dataset", values_from = "TFR") %>% 
  pivot_longer(cols = 2:10, names_to = "Dataset", values_to = "Genealogy") %>%
  mutate(Error = Genealogy - `Whole Simulation` , 
         Relative_Error = (Error/`Whole Simulation`)*100,
         Type = "DiM") %>% 
  select(-c(Genealogy,`Whole Simulation`)) 

# Mean of differences
MoD_TFR_Exp2 <- bind_rows(TFR_whole, TFR_dir_wod, TFR_anc_z, TFR_anc_au, TFR_anc_gau, TFR_anc_ggau, 
                          TFR_anc_gggau, TFR_anc_ggggau, TFR_anc_gggggau, TFR_anc_ggggggau) %>%
  filter(Year > 1750) %>% 
  # There can be rates of 0 and NaN values as the genealogist (most recent generation) are at least 18 years old. 
  # This happens after 2004
  filter(TFR != 0 & !is.nan(TFR)) %>% 
  pivot_wider(id_cols = c(Year, Sim_id), names_from = "Dataset", values_from = "TFR") %>% 
  pivot_longer(cols = 4:12, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation`, 
         Relative_Error = (Error/`Whole Simulation`)*100) %>% 
  group_by(Year, Dataset) %>% 
  summarise(Error = mean(Error, na.rm = T), 
          Relative_Error = mean(Relative_Error, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Type = "MoD")

# Bind both error measures and save the data frame
error_TFR_exp2 <- bind_rows(DiM_TFR_Exp2, MoD_TFR_Exp2)
save(error_TFR_exp2, file = "Measures/error_TFR_exp2.RData")

# Absolute Error  
error_TFR_exp2 %>% 
  ggplot(aes(x = Year, y = Error, group = Dataset, colour = Type)) +
  facet_wrap(. ~ Dataset)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp2_TFR_Error.jpeg", width=17, height=9, dpi=300)

# Relative Error
error_TFR_exp2 %>% 
  ggplot(aes(x = Year, y = Relative_Error, group = Dataset, colour = Type)) +
  facet_wrap(. ~ Dataset)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp2_TFR_Rel_Error.jpeg", width=17, height=9, dpi=300)

# Retrieve the year with the smallest error (close to 0)
error_TFR_exp2 %>%
  filter(Dataset == "8. + 6x-Great-Aunts/Uncles") %>% 
  filter(Error == max(Error))

# Check minimum and maximum values of bias in TFR before 1834
error_TFR_exp2 %>%
  filter(Year < 1834) %>% 
  filter(Dataset == "8. + 6x-Great-Aunts/Uncles") %>% 
  pull(Error) %>% # -0.17983396 -0.05247996
  range()

# Check mean value of bias and relative bias in TFR before 1834
error_TFR_exp2 %>% 
  filter(Year < 1851) %>% 
  filter(Dataset == "8. + 6x-Great-Aunts/Uncles") %>% 
  #pull(Error) %>% 
  pull(Relative_Error) %>% # -2.438331
  mean()

# Check minimum and maximum values of bias in TFR after 1834 # Not sure why using this number
error_TFR_exp2 %>%
  filter(Year > 1834) %>% 
  filter(Dataset == "8. + 6x-Great-Aunts/Uncles") %>% 
  pull(Error) %>% 
  range()

# Check maximum values of relative bias in TFR
error_TFR_exp2 %>%
  filter(Dataset == "8. + 6x-Great-Aunts/Uncles") %>% 
  filter(Relative_Error == min(Relative_Error))

# Check bias and relative bias in TFR in 2022 (year with the greatest error)
error_TFR_exp2 %>%
  filter(Year == 2022) %>% 
  filter(Dataset == "8. + 6x-Great-Aunts/Uncles")


# Life Expectancy at birth ----
# Estimate life expectancy at birth from asmr 1x1 for the different genealogical subsets ----

# Estimate age-specific mortality rates 1x1 from the subset of all direct ancestors and their offspring
asmr_anc_off_1 <- map_dfr(anc_off, ~ estimate_mortality_rates_mod(opop = .x,
                                                              final_sim_year = 2022, #[Jan-Dec]
                                                              year_min = 1750, # Closed
                                                              year_max = 2023, # Open )
                                                              year_group = 1,
                                                              age_max_mort = 110, # Open )
                                                              age_group = 1), # [,)
                          .id = "Sim_id") 
save(asmr_anc_off_1, file = "Measures/asmr_anc_off_1.RData")

# Compute life tables from the subset of all direct ancestors and their offspring
lt_anc_off <- lt_socsim_sims(asmr_socsim_sims = asmr_anc_off_1)
save(lt_anc_off, file = "Measures/lt_anc_off.RData")

# Estimate age-specific mortality rates 1x1 from the subset of direct ancestors, offspring above and siblings
asmr_anc_z_1 <- map_dfr(anc_z, ~ estimate_mortality_rates_mod(opop = .x,
                                                          final_sim_year = 2022, #[Jan-Dec]
                                                          year_min = 1750, # Closed
                                                          year_max = 2023, # Open )
                                                          year_group = 1,
                                                          age_max_mort = 110, # Open )
                                                          age_group = 1), # [,)
                        .id = "Sim_id") 
save(asmr_anc_z_1, file = "Measures/asmr_anc_z_1.RData")

# Compute life tables from the subset of direct ancestors, offspring above and siblings
lt_anc_z <- lt_socsim_sims(asmr_socsim_sims = asmr_anc_z_1)
save(lt_anc_z, file = "Measures/lt_anc_z.RData")

# Estimate age-specific mortality rates 1x1 from the subset of direct ancestors, offspring above and aunts/uncles
asmr_anc_au_1 <- map_dfr(anc_au, ~ estimate_mortality_rates_mod(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed
                                                            year_max = 2023, # Open )
                                                            year_group = 1,
                                                            age_max_mort = 110, # Open )
                                                            age_group = 1), # [,)
                         .id = "Sim_id") 
save(asmr_anc_au_1, file = "Measures/asmr_anc_au_1.RData")

# Compute life tables from the subset of direct ancestors, offspring above and aunts/uncles
lt_anc_au <- lt_socsim_sims(asmr_socsim_sims = asmr_anc_au_1)
save(lt_anc_au, file = "Measures/lt_anc_au.RData")

# Estimate age-specific mortality rates 1x1 from the subset of direct ancestors, offspring above and great-aunts/uncles
asmr_anc_gau_1 <- map_dfr(anc_gau, ~ estimate_mortality_rates_mod(opop = .x,
                                                              final_sim_year = 2022, #[Jan-Dec]
                                                              year_min = 1750, # Closed
                                                              year_max = 2023, # Open )
                                                              year_group = 1,
                                                              age_max_mort = 110, # Open )
                                                              age_group = 1), # [,)
                          .id = "Sim_id") 
save(asmr_anc_gau_1, file = "Measures/asmr_anc_gau_1.RData")

# Compute life tables from the subset of direct ancestors, offspring above and great-aunts/uncles, 
lt_anc_gau <- lt_socsim_sims(asmr_socsim_sims = asmr_anc_gau_1)
save(lt_anc_gau, file = "Measures/lt_anc_gau.RData")

# Estimate age-specific mortality rates 1x1 from the subset of direct ancestors, offspring above and 2x-great-aunts/uncles
asmr_anc_ggau_1 <- map_dfr(anc_ggau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                final_sim_year = 2022, #[Jan-Dec]
                                                                year_min = 1750, # Closed
                                                                year_max = 2023, # Open )
                                                                year_group = 1,
                                                                age_max_mort = 110, # Open )
                                                                age_group = 1), # [,)
                           .id = "Sim_id") 
save(asmr_anc_ggau_1, file = "Measures/asmr_anc_ggau_1.RData")

# Compute life tables from the subset of direct ancestors, offspring above and 2x-great-aunts/uncles, 
lt_anc_ggau <- lt_socsim_sims(asmr_socsim_sims = asmr_anc_ggau_1)
save(lt_anc_ggau, file = "Measures/lt_anc_ggau.RData")

# Estimate age-specific mortality rates 1x1 from the subset of direct ancestors, offspring above and 3x-great-aunts/uncles
asmr_anc_gggau_1 <- map_dfr(anc_gggau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                  final_sim_year = 2022, #[Jan-Dec]
                                                                  year_min = 1750, # Closed
                                                                  year_max = 2023, # Open )
                                                                  year_group = 1,
                                                                  age_max_mort = 110, # Open )
                                                                  age_group = 1), # [,)
                            .id = "Sim_id") 
save(asmr_anc_gggau_1, file = "Measures/asmr_anc_gggau_1.RData")

# Compute life tables from the subset of direct ancestors, offspring above and 3x-great-aunts/uncles, 
lt_anc_gggau <- lt_socsim_sims(asmr_socsim_sims = asmr_anc_gggau_1)
save(lt_anc_gggau, file = "Measures/lt_anc_gggau.RData")


# Estimate age-specific mortality rates 1x1 from the subset of direct ancestors, offspring above and 4x-great-aunts/uncles
asmr_anc_ggggau_1 <- map_dfr(anc_ggggau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                    final_sim_year = 2022, #[Jan-Dec]
                                                                    year_min = 1750, # Closed
                                                                    year_max = 2023, # Open )
                                                                    year_group = 1,
                                                                    age_max_mort = 110, # Open )
                                                                    age_group = 1), # [,)
                             .id = "Sim_id") 
save(asmr_anc_ggggau_1, file = "Measures/asmr_anc_ggggau_1.RData")

# Compute life tables from the subset of direct ancestors, offspring above and 4x-great-aunts/uncles, 
lt_anc_ggggau <- lt_socsim_sims(asmr_socsim_sims = asmr_anc_ggggau_1)
save(lt_anc_ggggau, file = "Measures/lt_anc_ggggau.RData")


# Estimate age-specific mortality rates 1x1 from the subset of direct ancestors, offspring above and 5x-great-aunts/uncles
asmr_anc_gggggau_1 <- map_dfr(anc_gggggau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                      final_sim_year = 2022, #[Jan-Dec]
                                                                      year_min = 1750, # Closed
                                                                      year_max = 2023, # Open )
                                                                      year_group = 1,
                                                                      age_max_mort = 110, # Open )
                                                                      age_group = 1), # [,)
                              .id = "Sim_id") 
save(asmr_anc_gggggau_1, file = "Measures/asmr_anc_gggggau_1.RData")

# Compute life tables from the subset of direct ancestors, offspring above and 5x-great-aunts/uncles, 
lt_anc_gggggau <- lt_socsim_sims(asmr_socsim_sims = asmr_anc_gggggau_1)
save(lt_anc_gggggau, file = "Measures/lt_anc_gggggau.RData")


# Estimate age-specific mortality rates 1x1 from the subset of direct ancestors, offspring above and 6x-great-aunts/uncles
asmr_anc_ggggggau_1 <- map_dfr(anc_ggggggau, ~ estimate_mortality_rates_mod(opop = .x,
                                                                        final_sim_year = 2022, #[Jan-Dec]
                                                                        year_min = 1750, # Closed
                                                                        year_max = 2023, # Open )
                                                                        year_group = 1,
                                                                        age_max_mort = 110, # Open )
                                                                        age_group = 1), # [,)
                               .id = "Sim_id") 
save(asmr_anc_ggggggau_1, file = "Measures/asmr_anc_ggggggau_1.RData")

# Compute life tables from the subset of direct ancestors, offspring above and 6x-great-aunts/uncles, 
lt_anc_ggggggau <- lt_socsim_sims(asmr_socsim_sims = asmr_anc_ggggggau_1)
save(lt_anc_ggggggau, file = "Measures/lt_anc_ggggggau.RData")

# Load and wrangle life tables for plotting ----

# Load life tables from each whole SOCSIM simulation
load("Measures/lt_10.RData")
# Load life tables from the subset of direct ancestors (without duplicates)
load("Measures/with_distinct/lt_dir_wod.RData")
# Load life tables from subset of direct ancestors and siblings
load("Measures/lt_anc_z.RData")
# Load life tables from subset of direct ancestors, offspring above and aunts/uncles
load("Measures/lt_anc_au.RData")
# Load life tables from subset of direct ancestors, offspring above and great-aunts/uncles
load("Measures/lt_anc_gau.RData")
# Load life tables from subset of direct ancestors, offspring above and 2x-great-aunts/uncles
load("Measures/lt_anc_ggau.RData")
# Load life tables from subset of direct ancestors, offspring above and 3x-great-aunts/uncles
load("Measures/lt_anc_gggau.RData")
# Load life tables from subset of direct ancestors, offspring above and 4x-great-aunts/uncles
load("Measures/lt_anc_ggggau.RData")
# Load life tables from subset of direct ancestors, offspring above and 5x-great-aunts/uncles
load("Measures/lt_anc_gggggau.RData")
# Load life tables from subset of direct ancestors, offspring above and 6x-great-aunts/uncles
load("Measures/lt_anc_ggggggau.RData")
# Load life tables from subset of all direct ancestors and their offspring
load("Measures/lt_anc_off.RData")

# Load asmr 1x1 for the 10 simulations, calculated on 2_Compare_Input_Output
load("Measures/asmr_10_1.RData")
# Year breaks. Extract all the unique numbers from the intervals 
year_breaks_mort_1 <- unique(as.numeric(str_extract_all(asmr_10_1$year, "\\d+", simplify = T)))

# Year range to filter data
year_range_mort_1 <- min(year_breaks_mort_1):max(year_breaks_mort_1-1)

# Whole SOCSIM simulation
lt_whole2 <-  lt_10 %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Whole Simulation",
         Rate = "e0")  %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# Direct ancestors without duplicates
lt_dir_wod2 <- lt_dir_wod %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct Ancestors",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# Direct ancestors and siblings
lt_anc_z2 <- lt_anc_z %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "1. + Siblings",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# Direct ancestors, offspring above and aunts/uncles
lt_anc_au2 <- lt_anc_au %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "2. + Aunts/Uncles",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# Direct ancestors, offspring above and great-aunts/uncles
lt_anc_gau2 <- lt_anc_gau %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "3. + Great-Aunts/Uncles",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# Direct ancestors, offspring above and 2x-great-aunts/uncles
lt_anc_ggau2 <- lt_anc_ggau %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "4. + 2x-Great-Aunts/Uncles",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# Direct ancestors, offspring above and 3x-great-aunts/uncles
lt_anc_gggau2 <- lt_anc_gggau %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "5. + 3x-Great-Aunts/Uncles",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# Direct ancestors, offspring above and 4x-great-aunts/uncles
lt_anc_ggggau2 <- lt_anc_ggggau %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "6. + 4x-Great-Aunts/Uncles",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# Direct ancestors, offspring above and 5x-great-aunts/uncles
lt_anc_gggggau2 <- lt_anc_gggggau %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "7. + 5x-Great-Aunts/Uncles",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# Direct ancestors, offspring above and 6x-great-aunts/uncles
lt_anc_ggggggau2 <- lt_anc_ggggggau %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "8. + 6x-Great-Aunts/Uncles",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors, offspring above and their offspring
lt_anc_off2 <- lt_anc_off %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct Ancestors and their Offspring",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)


# Plot the estimates of life expectancy at birth
bind_rows(lt_whole2, lt_dir_wod2, lt_anc_z2, lt_anc_au2, lt_anc_gau2, lt_anc_ggau2, 
          lt_anc_gggau2, lt_anc_ggggau2, lt_anc_gggggau2, lt_anc_ggggggau2)  %>%
  filter(Age == 0) %>%
  mutate(Sex = ifelse(sex == "female", "Female", "Male"),
         Dataset = factor(Dataset, levels =  c("Direct Ancestors", "1. + Siblings", "2. + Aunts/Uncles", "3. + Great-Aunts/Uncles", "4. + 2x-Great-Aunts/Uncles", 
                                               "5. + 3x-Great-Aunts/Uncles", "6. + 4x-Great-Aunts/Uncles", "7. + 5x-Great-Aunts/Uncles", "8. + 6x-Great-Aunts/Uncles",
                                               "Whole Simulation"))) %>%
  group_by(Year, Sex, Dataset) %>% 
  summarise(ex = mean(ex, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year, y = ex, group = Dataset))+
  geom_line(aes(colour = Dataset), linewidth = 1.3)+
  scale_color_manual(values = c("#7A7500", "#FBD724FF", "#FEB82CFF", "#F1804DFF", "#E4695EFF", 
                                "#D45270FF", "#C23C81FF", "#AC2694FF", "#75007A", "#007A75"))+
  facet_wrap(~Sex) +
  theme_graphs() +
  labs(y = "Life expectancy at birth")
ggsave(file="Graphs/Socsim_Exp2_e0.jpeg", width=17, height=9, dpi=300)

#----------------------------------------------------------------------------------------------------
# Figure for EPC presentation

# Plot the estimates of life expectancy at birth
bind_rows(lt_whole2, lt_dir_wod2, lt_anc_off2)  %>%
  filter(Age == 0) %>%
  mutate(Sex = ifelse(sex == "female", "Female", "Male")) %>%
  group_by(Year, Sex, Dataset) %>% 
  summarise(ex = mean(ex, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year, y = ex, group = Dataset))+
  geom_line(aes(colour = Dataset), linewidth = 1.3)+
  scale_color_manual(values = c("#7A7500", "#75007A", "#007A75"))+
  facet_wrap(~Sex) +
  theme_graphs() +
  labs(y = "Life expectancy at birth")
ggsave(file="Graphs/Socsim_Exp2_e0_grp.jpeg", width=17, height=9, dpi=300)
#----------------------------------------------------------------------------------------------------
# Summary measure of error in e0 ----

# Difference in means
DiM_e0_Exp2 <- bind_rows(lt_whole2, lt_dir_wod2, lt_anc_z2, lt_anc_au2, lt_anc_gau2, lt_anc_ggau2, 
                         lt_anc_gggau2, lt_anc_ggggau2, lt_anc_gggggau2, lt_anc_ggggggau2) %>%
  filter(Year > 1750 & Age == 0) %>% 
  group_by(Year, sex, Dataset) %>% 
  summarise(ex = mean(ex, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(Year:sex), names_from = "Dataset", values_from = "ex") %>% 
  pivot_longer(cols = 3:11, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation` , 
         Relative_Error = (Error/`Whole Simulation`)*100,
         Type = "DiM") %>% 
  select(-c(Genealogy,`Whole Simulation`)) 

# Mean of differences
MoD_e0_Exp2 <- bind_rows(lt_whole2, lt_dir_wod2, lt_anc_z2, lt_anc_au2, lt_anc_gau2, lt_anc_ggau2, 
                         lt_anc_gggau2, lt_anc_ggggau2, lt_anc_gggggau2, lt_anc_ggggggau2) %>%
  filter(Year > 1750 & Age == 0) %>% 
  pivot_wider(id_cols = c(Year, sex, Sim_id), names_from = "Dataset", values_from = "ex") %>% 
  pivot_longer(cols = 5:13, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation`, 
         Relative_Error = (Error/`Whole Simulation`)*100) %>% 
  group_by(Year, sex, Dataset) %>% 
  summarise(Error = mean(Error, na.rm = T), 
          Relative_Error = mean(Relative_Error, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Type = "MoD")

# Bind both error measures and save the data frame
error_e0_exp2 <- bind_rows(DiM_e0_Exp2, MoD_e0_Exp2) 
save(error_e0_exp2, file = "Measures/error_e0_exp2.RData")

# Absolute error
error_e0_exp2 %>% 
  ggplot(aes(x = Year, y = Error, group = Dataset, colour = Type)) +
  facet_grid(Dataset ~ sex)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp2_e0_Error.jpeg", width=17, height=9, dpi=300)

# Relative error
error_e0_exp2 %>% 
  ggplot(aes(x = Year, y = Relative_Error, group = Dataset, colour = Type)) +
  facet_grid(sex ~ Dataset)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp2_e0_Rel_Error.jpeg", width=17, height=9, dpi=300)

# Check minimum and maximum values of bias in e0 
error_e0_exp2 %>% 
  filter(sex == "female" & Dataset == "8. + 6x-Great-Aunts/Uncles") %>% 
  pull(Error) %>%
  range()

## Check years when error is negative (so, underestimation) 
error_e0_exp2 %>% 
  filter(sex == "female" & Dataset == "8. + 6x-Great-Aunts/Uncles") %>% 
  filter(Error < 0) %>% 
  pull(Year)

# Check mean values of bias in e0 
error_e0_exp2 %>% 
  filter(sex == "female" & Dataset == "8. + 6x-Great-Aunts/Uncles") %>% 
  pull(Error) %>%
  mean()

# Check mean values of relative bias in e0 
error_e0_exp2 %>% 
  filter(sex == "female" & Dataset == "8. + 6x-Great-Aunts/Uncles") %>% 
  pull(Relative_Error) %>%
  mean()

#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

yrs_plot2 <- c(1750, 1800, 1850, 1900, 1950, 2000)

# Define the same y breaks for all plots
y_breaks_TFR <- c(0:5)
y_breaks_e0 <- c(20, 40, 60, 80)

## TFR and e0 (for females) from whole SOCSIM simulation and subsets of direct ancestors without and with their offspring
Summary_Exp2 <- 
bind_rows(TFR_whole %>% rename(Estimate = TFR),
          # There can be TFR of 0 and NaN values as the genealogist (most recent generation) are at least 18 years old. 
          TFR_dir_wod %>% rename(Estimate = TFR) %>% filter(Estimate != 0 & !is.nan(Estimate)),
          TFR_anc_off %>% rename(Estimate = TFR)) %>%  
  mutate(sex = "female") %>%  
  bind_rows(lt_whole2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_dir_wod2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_off2 %>% rename(Estimate = ex) %>% filter(Age == 0)) %>%
  filter(sex == "female") %>%
  group_by(Year, Dataset, Rate, sex) %>% 
  summarise(Estimate = mean(Estimate, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life Expectancy at Birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life Expectancy at Birth")), 
         Dataset = factor(Dataset, levels =  c("Direct Ancestors", "Direct Ancestors and their Offspring", "Whole Simulation"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), 
             aes(shape = Dataset), size = 11) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#7A7500", "#75007A", "#007A75"))+
  scale_shape_manual(values = c(19, 18, 46)) +
  scale_x_continuous(breaks = yrs_plot2)+
  facetted_pos_scales(y = list("Total Fertility Rate" = scale_y_continuous(breaks = y_breaks_TFR, 
                                                                           limits = c(0, NA)),
                               "Life Expectancy at Birth" =  scale_y_continuous(breaks = y_breaks_e0)))+
  theme_graphs() +
  theme(legend.justification = "left",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

# Save the plot
Summary_Exp2
ggsave(file="Graphs/Socsim_Exp2_TFR_e0.jpeg", width=17, height=9, dpi=300)

#----------------------------------------------------------------------------------------------------
## Plot combining age-specific rates and summary measures -----

plot_labs1 <- data.frame(Rate = c("Age-Specific Fertility Rates", "Age-Specific Mortality Rates"),
                         x = c(1,2),
                         y = c(0.23, 0.55),
                         labels = c("a","b"))
plot_labs2 <- data.frame(Rate = as.factor(c("Total Fertility Rate", "Life Expectancy at Birth")),
                         x = c(1755, 1755),
                         y = c(5.4, 84),
                         labels = c("c","d"))

By_Age_Exp2 + 
  geom_text(data = plot_labs1, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") + 
  theme(plot.margin = margin(0,0,1,0, "cm")) +
  Summary_Exp2 + 
  geom_text(data = plot_labs2, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") +
  plot_layout(ncol = 1)
ggsave(file="Final_Graphs/Final_Socsim_Exp2_Combined.jpeg", width=18, height=21, dpi=300)

#----------------------------------------------------------------------------------------------------
## Plots adding kin progressively For appendix ----

## Plotting ASFR and ASMR (for females) from whole SOCSIM simulation and subsets of direct ancestors and their offspring
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2 %>% rename(Estimate = ASFR),
          asfr_dir_wod2 %>% rename(Estimate = ASFR),
          asfr_anc_z2 %>% rename(Estimate = ASFR),
          asfr_anc_au2 %>% rename(Estimate = ASFR),
          asfr_anc_gau2 %>% rename(Estimate = ASFR),
          asfr_anc_ggau2 %>% rename(Estimate = ASFR),
          asfr_anc_gggau2 %>% rename(Estimate = ASFR),
          asfr_anc_ggggau2 %>% rename(Estimate = ASFR),
          asfr_anc_gggggau2 %>% rename(Estimate = ASFR),
          asfr_anc_ggggggau2 %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx),
            asmr_dir_wod2 %>% rename(Estimate = mx),
            asmr_anc_z2 %>% rename(Estimate = mx),
            asmr_anc_au2 %>% rename(Estimate = mx),
            asmr_anc_gau2 %>% rename(Estimate = mx),
            asmr_anc_ggau2 %>% rename(Estimate = mx),
            asmr_anc_gggau2 %>% rename(Estimate = mx),
            asmr_anc_ggggau2 %>% rename(Estimate = mx),
            asmr_anc_gggggau2 %>% rename(Estimate = mx),
            asmr_anc_ggggggau2 %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels),
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates"), 
         Dataset = ifelse(Dataset == "Direct Ancestors", "0. Direct Ancestors", Dataset),
         Dataset = factor(Dataset, levels =  c("0. Direct Ancestors", "1. + Siblings", "2. + Aunts/Uncles", "3. + Great-Aunts/Uncles", "4. + 2x-Great-Aunts/Uncles", 
                                               "5. + 3x-Great-Aunts/Uncles", "6. + 4x-Great-Aunts/Uncles", "7. + 5x-Great-Aunts/Uncles", "8. + 6x-Great-Aunts/Uncles",
                                               "Whole Simulation"))) %>%
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Dataset))+
  facet_wrap(Year ~ Rate, nrow = 3, ncol = 2, scales = "free") + 
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(size = 9)+ 
  scale_color_manual(values = c("#7A7500",  "#FBD724FF", "#FEB82CFF", "#F1804DFF", "#E4695EFF", 
                                "#D45270FF", "#C23C81FF", "#AC2694FF", "#75007A", "#007A75"))+
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10"),
                               ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10"),
                               ASFR = scale_y_continuous(breaks = y_breaks_asfr, limits = c(0, 0.2)),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10"))) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(x = "Age")
ggsave(file="Final_Graphs/App_Socsim_Exp2_ASFR_ASMR.jpeg", width=18, height=25, dpi=300)

## TFR and e0 (for females) from whole SOCSIM simulation and subsets of direct ancestors without and with their offspring
bind_rows(TFR_whole %>% rename(Estimate = TFR),
          # There can be TFR of 0 and NaN values as the genealogist (most recent generation) are at least 18 years old. 
          TFR_dir_wod %>% rename(Estimate = TFR) %>% filter(Estimate != 0 & !is.nan(Estimate)),
          TFR_anc_z %>% rename(Estimate = TFR) %>% filter(Estimate != 0 & !is.nan(Estimate)),
          TFR_anc_au %>% rename(Estimate = TFR) %>% filter(Estimate != 0 & !is.nan(Estimate)),
          TFR_anc_gau %>% rename(Estimate = TFR),
          TFR_anc_ggau %>% rename(Estimate = TFR),
          TFR_anc_gggau %>% rename(Estimate = TFR),
          TFR_anc_ggggau %>% rename(Estimate = TFR),
          TFR_anc_gggggau %>% rename(Estimate = TFR),
          TFR_anc_ggggggau %>% rename(Estimate = TFR)) %>%  
  mutate(sex = "female") %>%  
  bind_rows(lt_whole2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_dir_wod2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_z2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_au2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_gau2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_ggau2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_gggau2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_ggggau2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_gggggau2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_ggggggau2 %>% rename(Estimate = ex) %>% filter(Age == 0)) %>%
  filter(sex == "female") %>%
  group_by(Year, Dataset, Rate, sex) %>% 
  summarise(Estimate = mean(Estimate, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life Expectancy at Birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life Expectancy at Birth")), 
         Dataset = ifelse(Dataset == "Direct Ancestors", "0. Direct Ancestors", Dataset),
         Dataset = factor(Dataset, levels =  c("0. Direct Ancestors", "1. + Siblings", "2. + Aunts/Uncles", "3. + Great-Aunts/Uncles", "4. + 2x-Great-Aunts/Uncles", 
                                               "5. + 3x-Great-Aunts/Uncles", "6. + 4x-Great-Aunts/Uncles", "7. + 5x-Great-Aunts/Uncles", "8. + 6x-Great-Aunts/Uncles",
                                               "Whole Simulation"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#7A7500", "#FBD724FF", "#FEB82CFF", "#F1804DFF", "#E4695EFF", 
                                "#D45270FF", "#C23C81FF", "#AC2694FF", "#75007A", "#007A75"))+
  facetted_pos_scales(y = list("Total Fertility Rate" = scale_y_continuous(breaks = y_breaks_TFR, 
                                                                           limits = c(0, NA)),
                               "Life Expectancy at Birth" =  scale_y_continuous(breaks = y_breaks_e0)))+
  theme_graphs()
# Save the plot
ggsave(file="Final_Graphs/App_Socsim_Exp2_TFR_e0.jpeg", width=17, height=11, dpi=300)