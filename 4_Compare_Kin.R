#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Sweden - Trace and compare genealogical subsets of kin of a given ego(s) 
# U:/SOCSIM/SOCSIM_Genealogies/4_Compare_Kin.R

## Trace relevant ascending and lateral kin up to the 4th generation
# of a given ego(s) from a SOCSIM microsimulation for Sweden (1751-2022)
# and compare demographic measures from the whole simulation and the genealogical subsets

# Created by Liliana Calderon on 13-04-2022
# Last modified by Liliana Calderon on 10-07-2023

## NB: To run this code, it is necessary to have already run the scripts 
# 1_Run_Simulations.R and 3_Compare_Ancestors.R

#------------------------------------------------------------------------------------------------------
## General settings and functions ----

# Prevent scientific notation (useful for the rate calculation)
options(scipen=999999)

## Load packages 
library(tidyverse)
library(ggh4x)  # To facet scales-
library(svglite) # To save svg files
library(viridis)
library(rsocsim) # Functions to estimate rates

## Load theme for the graphs and to convert SOCSIM time
source("Functions_Graphs.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with asmr calculated with rsocsim::estimate_mortality_rates()
source("Functions_Life_Table.R")

## Load function to re-code the kin_type for the ancestors
source("Functions_Ancestors.R")

#------------------------------------------------------------------------------------------------------
## Trace kin (up to 4th degree) of people alive in 2023 ----

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
                             ~ rsocsim::retrieve_kin(opop = ..1, omar = ..2, pid = ..3,
                                                     KidsOf = KidsOf,
                                                     extra_kintypes = c("unclesaunts", "gunclesaunts", "firstcousins"),
                                                     kin_by_sex = F),
                             .id = "Sim_id")
# Time difference of 26.02852 mins for 10 simulations with initial population of 5000, with hetfert0

# Select relevant variables from each opop and add Sim_id to allow merging with the kin_10 df 
sims_opop_id <- map_dfr(sims_opop, ~ select(.x, c(pid, fem, dob, dod, mom, lborn, marid, mstat)), 
                        .id = "Sim_id")

# Format the kin data,  put into long format and add columns from each opop
# NB: The function in the package still have some duplicates that need to be removed
# With the kin_by_sex = F option, we get nieces and nephews in addition to niblings
kin_10 <- kin_10_lc %>% 
  mutate(ego_id = unlist(egos_samp_10)) %>% 
  pivot_longer(-c(Sim_id, ego_id), names_to = "kin_type", values_to = "pid") %>% 
  unnest_longer(kin_type:pid) %>%
  filter(!is.na(pid)) %>% 
  # Sometimes, a pid is included as more than one kin type, 
  # e.g. ggparents and gunclesaunts or siblings and firstcousins of the same ego
  ## because parents were included as (gg)parents and (g)unclesaunts
  distinct(Sim_id, ego_id, pid, .keep_all = TRUE) %>% 
  left_join(sims_opop_id, by = c("Sim_id", "pid"))

# Save the data frame with the relatives of 10 simulations samples
save(kin_10, file = "Subsets/kin_10.RData")

#------------------------------------------------------------------------------------------------------
## Merge direct ancestors with lateral kin of sample of egos alive in 2023 ----

# Load the data frame with the relatives of 10 simulations samples
load("Subsets/kin_10.RData")

# Load the data frame with the ancestors of the 10 simulations samples of egos alive in 2023
load("Subsets/ancestors_10.RData")

## Re-code the kin_type for the ancestors and remove duplicates (from the different family trees)
ancestors_10 <- ancestors_10 %>% 
  recode_ancestors() %>%
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all = TRUE) %>% 
  ungroup() 

# Merge the two data frames and remove kin types included in both (i.e.g, parents, gparents, ggparents)
# NB: This data still contains duplicates within each simulation, as someone can be a relative of different egos
anc_kin_10 <- bind_rows(ancestors_10, kin_10) %>% 
  group_by(Sim_id, ego_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup()

# Save the data frame with the ancestors and relatives of 10 simulations samples
save(anc_kin_10, file = "Subsets/anc_kin_10.RData")

#----------------------------------------------------------------------------------------------------
## Age-Specific Fertility and Mortality rates, 5x5  -----
# Retrieve and compare rates derived from the whole simulation with those from the genealogical subsets
# with direct and lateral kin up to the 4th degree of consanguinity

load("Subsets/anc_kin_10.RData")

## This repeated code could be improved using a function to calculate the rates for each subset, 
# c.f. trial on "progress/Functions_Rates_Kin.R". 
# The function worked but I did not manage to assign a different name to each data frame

##  Estimate ASFR and ASMR for the genealogical subsets with direct ancestors and collateral kin

## Direct ancestors and siblings ----
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

# Estimate age-specific fertility rates for the subset of direct ancestors and siblings
asfr_anc_z <- map_dfr(anc_z, ~ estimate_fertility_rates(opop = .x,
                                                        final_sim_year = 2022, #[Jan-Dec]
                                                        year_min = 1750, # Closed [
                                                        year_max = 2020, # Open )
                                                        year_group = 5, 
                                                        age_min_fert = 10, # Closed [
                                                        age_max_fert = 55, # Open )
                                                        age_group = 5), # [,)
                      .id = "Sim_id") 
save(asfr_anc_z, file = "Measures/asfr_anc_z.RData")

# Estimate age-specific mortality rates for the subset of direct ancestors and siblings
asmr_anc_z <- map_dfr(anc_z, ~ estimate_mortality_rates(opop = .x,
                                                        final_sim_year = 2022, #[Jan-Dec]
                                                        year_min = 1750, # Closed
                                                        year_max = 2020, # Open )
                                                        year_group = 5,
                                                        age_max_mort = 110, # Open )
                                                        age_group = 5), # [,)
                      .id = "Sim_id") 
save(asmr_anc_z, file = "Measures/asmr_anc_z.RData")


## Direct ancestors and aunts/uncles 8% ----
type_anc_au <- c("ego", "parents", "gparents", "ggparents", 
                 "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                 "unclesaunts")

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_au <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_au) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)

# Estimate age-specific fertility rates for the subset of direct ancestors and aunts/uncles
asfr_anc_au <- map_dfr(anc_au, ~ estimate_fertility_rates(opop = .x,
                                                          final_sim_year = 2022, #[Jan-Dec]
                                                          year_min = 1750, # Closed [
                                                          year_max = 2020, # Open )
                                                          year_group = 5, 
                                                          age_min_fert = 10, # Closed [
                                                          age_max_fert = 55, # Open )
                                                          age_group = 5), # [,)
                       .id = "Sim_id") 
save(asfr_anc_au, file = "Measures/asfr_anc_au.RData")

# Estimate age-specific mortality rates for the subset of direct ancestors and aunts/uncles
asmr_anc_au <- map_dfr(anc_au, ~ estimate_mortality_rates(opop = .x,
                                                          final_sim_year = 2022, #[Jan-Dec]
                                                          year_min = 1750, # Closed
                                                          year_max = 2020, # Open )
                                                          year_group = 5,
                                                          age_max_mort = 110, # Open )
                                                          age_group = 5), # [,)
                       .id = "Sim_id") 
save(asmr_anc_au, file = "Measures/asmr_anc_au.RData")


## Direct ancestors and first cousins, 14% ----
type_anc_k <- c("ego", "parents", "gparents", "ggparents", 
                "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                "firstcousins") 

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_k <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_k) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)

# Estimate age-specific fertility rates for the subset of direct ancestors and cousins
asfr_anc_k <- map_dfr(anc_k, ~ estimate_fertility_rates(opop = .x,
                                                        final_sim_year = 2022, #[Jan-Dec]
                                                        year_min = 1750, # Closed [
                                                        year_max = 2020, # Open )
                                                        year_group = 5, 
                                                        age_min_fert = 10, # Closed [
                                                        age_max_fert = 55, # Open )
                                                        age_group = 5), # [,)
                      .id = "Sim_id") 
save(asfr_anc_k, file = "Measures/asfr_anc_k.RData")

# Estimate age-specific mortality rates for the subset of direct ancestors and cousins
asmr_anc_k <- map_dfr(anc_k, ~ estimate_mortality_rates(opop = .x,
                                                        final_sim_year = 2022, #[Jan-Dec]
                                                        year_min = 1750, # Closed
                                                        year_max = 2020, # Open )
                                                        year_group = 5,
                                                        age_max_mort = 110, # Open )
                                                        age_group = 5), # [,)
                      .id = "Sim_id") 
save(asmr_anc_k, file = "Measures/asmr_anc_k.RData")


## Direct ancestors and great-aunt/uncles 22% ----
type_anc_gau <- c("ego", "parents", "gparents", "ggparents", 
                  "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                  "gunclesaunts") 

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_gau <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_gau) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)

# Estimate age-specific fertility rates for the subset of direct ancestors and great-aunt/uncles
asfr_anc_gau <- map_dfr(anc_gau, ~ estimate_fertility_rates(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed [
                                                            year_max = 2020, # Open )
                                                            year_group = 5, 
                                                            age_min_fert = 10, # Closed [
                                                            age_max_fert = 55, # Open )
                                                            age_group = 5), # [,)
                        .id = "Sim_id") 
save(asfr_anc_gau, file = "Measures/asfr_anc_gau.RData")

# Estimate age-specific mortality rates for the subset of direct ancestors and great-aunt/uncles
asmr_anc_gau <- map_dfr(anc_gau, ~ estimate_mortality_rates(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed
                                                            year_max = 2020, # Open )
                                                            year_group = 5,
                                                            age_max_mort = 110, # Open )
                                                            age_group = 5), # [,)
                        .id = "Sim_id") 
save(asmr_anc_gau, file = "Measures/asmr_anc_gau.RData")

## Direct ancestors and spouse 1% and children 3%  ----

type_anc_sc <- c("ego", "parents", "gparents", "ggparents", 
                 "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                 "spouse", "children") 

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_sc <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_sc) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)

# Estimate age-specific fertility rates for the subset of direct ancestors and spouse, children
asfr_anc_sc <- map_dfr(anc_sc, ~ estimate_fertility_rates(opop = .x,
                                                          final_sim_year = 2022, #[Jan-Dec]
                                                          year_min = 1750, # Closed [
                                                          year_max = 2020, # Open )
                                                          year_group = 5, 
                                                          age_min_fert = 10, # Closed [
                                                          age_max_fert = 55, # Open )
                                                          age_group = 5), # [,)
                       .id = "Sim_id") 
save(asfr_anc_sc, file = "Measures/asfr_anc_sc.RData")

# Estimate age-specific mortality rates for the subset of direct ancestors and spouse, children
asmr_anc_sc <- map_dfr(anc_sc, ~ estimate_mortality_rates(opop = .x,
                                                          final_sim_year = 2022, #[Jan-Dec]
                                                          year_min = 1750, # Closed
                                                          year_max = 2020, # Open )
                                                          year_group = 5,
                                                          age_max_mort = 110, # Open )
                                                          age_group = 5), # [,)
                       .id = "Sim_id") 
save(asmr_anc_sc, file = "Measures/asmr_anc_sc.RData")


## Direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles, spouse and children  ----

type_anc_zaukgausc <- c("ego", "parents", "gparents", "ggparents", 
                        "gggparents", "ggggparents", "gggggparents", "ggggggparents", "gggggggparents", 
                        "siblings", "unclesaunts", "firstcousins", "gunclesaunts", 
                        "spouse", "children") 

# Create a list of data frames with opop of selected kin keeping only unique pids for each simulation
anc_zaukgausc <- anc_kin_10 %>% 
  filter(kin_type %in% type_anc_zaukgausc) %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)

# Estimate age-specific fertility rates for the subset of direct ancestors, siblings, aunts/uncles, cousins, great-aunt/uncles, 
# spouse, children
asfr_anc_zaukgausc <- map_dfr(anc_zaukgausc, ~ estimate_fertility_rates(opop = .x,
                                                                        final_sim_year = 2022, #[Jan-Dec]
                                                                        year_min = 1750, # Closed [
                                                                        year_max = 2020, # Open )
                                                                        year_group = 5, 
                                                                        age_min_fert = 10, # Closed [
                                                                        age_max_fert = 55, # Open )
                                                                        age_group = 5), # [,)
                              .id = "Sim_id") 
save(asfr_anc_zaukgausc, file = "Measures/asfr_anc_zaukgausc.RData")

# Estimate age-specific mortality rates for the subset of direct ancestors, siblings, aunts/uncles, cousins, great-aunt/uncles, 
# spouse, children
asmr_anc_zaukgausc <- map_dfr(anc_zaukgausc, ~ estimate_mortality_rates(opop = .x,
                                                                        final_sim_year = 2022, #[Jan-Dec]
                                                                        year_min = 1750, # Closed
                                                                        year_max = 2020, # Open )
                                                                        year_group = 5,
                                                                        age_max_mort = 110, # Open )
                                                                        age_group = 5), # [,)
                              .id = "Sim_id") 
save(asmr_anc_zaukgausc, file = "Measures/asmr_anc_zaukgausc.RData")


#----------------------------------------------------------------------------------------------------
## Plot estimates from genealogical subsets of direct and lateral kin----

# Load ASFR and ASMR for the subset of "direct" family trees without duplicates
load("Measures/asfr_dir_wod.RData")
load("Measures/asmr_dir_wod.RData")

# Load asfr for the subset of direct ancestors and siblings
load("Measures/asfr_anc_z.RData")
# Load asmr for the subset of direct ancestors and siblings
load("Measures/asmr_anc_z.RData")

# Load asfr for the subset of direct ancestors and aunts/uncles
load("Measures/asfr_anc_au.RData")
# Load asmr for the subset of  direct ancestors and aunts/uncles
load("Measures/asmr_anc_au.RData")

# Load asfr for the subset of direct ancestors and cousins
load("Measures/asfr_anc_k.RData")
# Load asmr for the subset of direct ancestors and cousins
load("Measures/asmr_anc_k.RData")

# Load asfr for the subset of direct ancestors and great-aunts/uncles
load("Measures/asfr_anc_gau.RData")
# Load asmr for the subset of direct ancestors and great-aunts/uncles
load("Measures/asmr_anc_gau.RData")

# Load asfr for the subset of direct ancestors, spouse and children
load("Measures/asfr_anc_sc.RData")
# Load asmr for the subset of direct ancestors, spouse and children
load("Measures/asmr_anc_sc.RData")

# Load asfr for the subset of direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles
# spouse and children
load("Measures/asfr_anc_zaukgausc.RData")
# Load asmr for the subset of direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles
# spouse and children
load("Measures/asmr_anc_zaukgausc.RData")

# Choose years to plot (in intervals).
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_dir_wod$age)

# Function to plot asfr and asmr from each subset, for women
plot_asfr_asmr <- function(asfr, asmr, yrs_plot = yrs_plot, age_levels = age_levels) {
  bind_rows(asfr %>% 
              mutate(rate = "ASFR",
                     sex = "female"),
            asmr %>% 
              mutate(rate = "ASMR") %>% 
              filter(sex == "female")) %>% 
    # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
    filter(socsim !=0 & !is.infinite(socsim) & !is.nan(socsim)) %>% 
    filter(year %in% yrs_plot) %>% 
    mutate(age = factor(as.character(age), levels = age_levels)) %>% 
    ggplot(aes(x = age, y = socsim, group = interaction(Sim_id, year), colour = year)) +
    geom_line(linewidth = 1) +
    facet_wrap(. ~ rate, scales = "free") + 
    theme_graphs()  +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                                 ASMR =  scale_y_continuous(trans = "log10"))) +
    scale_color_viridis(option = "D", discrete = T, direction = -1)+
    labs(x = "Age", y = "Estimate")
}

## Plot ASFR and ASMR (for women) for the subset of "direct" family trees
plot_asfr_asmr(asfr_dir_wod, asmr_dir_wod, yrs_plot, age_levels)

# Plot ASFR and ASMR (for women) for the subset of direct ancestors and siblings
plot_asfr_asmr(asfr_anc_z, asmr_anc_z, yrs_plot, age_levels)

# Plot ASFR and ASMR (for women) for the subset of direct ancestors and aunts/uncles
plot_asfr_asmr(asfr_anc_au, asmr_anc_au, yrs_plot, age_levels)

# Plot ASFR and ASMR (for women) for the subset of direct ancestors and cousins
plot_asfr_asmr(asfr_anc_k, asmr_anc_k, yrs_plot, age_levels)

# Plot ASFR and ASMR (for women) for the subset of direct ancestors and great-aunts/uncles
plot_asfr_asmr(asfr_anc_gau, asmr_anc_gau, yrs_plot, age_levels)

# Plot ASFR and ASMR (for women) for the subset of direct ancestors and spouse and children
plot_asfr_asmr(asfr_anc_sc, asmr_anc_sc, yrs_plot, age_levels)

# Plot ASFR and ASMR (for women) for the subset of direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles
# spouse and children
plot_asfr_asmr(asfr_anc_zaukgausc, asmr_anc_zaukgausc, yrs_plot, age_levels)

#----------------------------------------------------------------------------------------------------
## Comparison of whole simulation with subsets of direct and extended kin ----

# ASFR ----

# Load mean ASFR 5x5 from the 10 simulations, calculated on 2_Compare_Input_Output
load("Measures/asfr_whole.RData")

## Calculate the mean of the different simulations and add relevant columns

# Whole SOCSIM simulation
asfr_whole2 <- asfr_whole %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Whole Simulation", 
         Rate = "ASFR") 

# "Direct" family trees without duplicates
asfr_dir_wod2 <- asfr_dir_wod %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct Ancestors (DA)",
         Rate = "ASFR") 

# Direct ancestors and siblings
asfr_anc_z2 <- asfr_anc_z %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "DA + Siblings",
         Rate = "ASFR") 

# Direct ancestors, siblings and aunts/uncles
asfr_anc_au2 <- asfr_anc_au %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "DA + Aunts/Uncles",
         Rate = "ASFR") 

# Direct ancestors, siblings, aunts/uncles and cousins
asfr_anc_k2 <- asfr_anc_k %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "DA + Cousins",
         Rate = "ASFR") 

# Direct ancestors, siblings, aunts/uncles, cousins, and great-aunts/uncles
asfr_anc_gau2 <- asfr_anc_gau %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "DA + Great-Aunts/Uncles",
         Rate = "ASFR") 

# Direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles
# spouse and children
asfr_anc_sc2 <- asfr_anc_sc %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "DA + Spouse and Children",
         Rate = "ASFR") 

# Direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles
# spouse and children
asfr_anc_zaukgausc2 <- asfr_anc_zaukgausc %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "ASFR") 


## Plot ASFR from whole SOCSIM simulation and subsets of direct and extended family trees without duplicates

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2, asfr_dir_wod2, asfr_anc_z2, asfr_anc_z2, asfr_anc_au2, 
          asfr_anc_k2, asfr_anc_gau2, asfr_anc_sc2)  %>% 
  filter(year %in% yrs_plot & !is.nan(ASFR)) %>% 
  mutate(Dataset = factor(Dataset, levels =  c("Direct Ancestors (DA)", "DA + Siblings", "DA + Aunts/Uncles", "DA + Cousins", "DA + Great-Aunts/Uncles", 
                                               "DA + Spouse and Children", "Whole Simulation"))) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset)))+
  facet_grid(year ~ . , scales = "free") +
  geom_line(aes(colour = year), linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(colour = year, shape = Dataset), size = 11)+ 
  scale_color_manual(values = c("#B72779", "#2779B7"))+
  scale_shape_manual(values = c(19, 10, 8, 2, 11, 7, 46)) +
  theme_graphs()
# labs(title = "Age-Specific Fertility Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and subsets of "direct" and "extended" family trees") 
ggsave(file="Graphs/Socsim_Exp2_ASFR.jpeg", width=17, height=14, dpi=200)

# ASMR ----

# Load mean ASMR rates 5x5 from the 10 simulations, calculated on 2_Compare_Input_Output
load("Measures/asmr_whole.RData")

# Whole SOCSIM simulation
asmr_whole2 <- asmr_whole %>% 
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
         Dataset = "Direct Ancestors (DA)",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Direct ancestors and siblings
asmr_anc_z2 <- asmr_anc_z %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "DA + Siblings",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Direct ancestors and aunts/uncles
asmr_anc_au2 <- asmr_anc_au %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "DA + Aunts/Uncles",
         Rate = "ASMR") %>%   
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Direct ancestors and cousins
asmr_anc_k2 <- asmr_anc_k %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "DA + Cousins",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Direct ancestors and great-aunts/uncles
asmr_anc_gau2 <- asmr_anc_gau %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "DA + Great-Aunts/Uncles",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Direct ancestors,spouse and children
asmr_anc_sc2 <- asmr_anc_sc %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "DA + Spouse and Children",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles, spouse and children
asmr_anc_zaukgausc2 <- asmr_anc_zaukgausc %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)


## Plotting ASMR from whole SOCSIM simulation and subsets of "direct" and "extended" family trees without duplicates

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1900,1905)", "[2000,2005)") 

bind_rows(asmr_whole2, asmr_dir_wod2, asmr_anc_z2, asmr_anc_z2, asmr_anc_au2, 
          asmr_anc_k2, asmr_anc_gau2, asmr_anc_sc2)  %>% 
  rename(Year = year) %>% 
  filter(Year %in% yrs_plot) %>%  
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>% 
  mutate(Dataset = factor(Dataset, levels =  c("Direct Ancestors (DA)", "DA + Siblings", "DA + Aunts/Uncles", "DA + Cousins", "DA + Great-Aunts/Uncles", 
                                               "DA + Spouse and Children",   "Whole Simulation"))) %>% 
  ggplot(aes(x = age, y = mx, group = interaction(Year, Dataset)))+
  facet_grid(Year ~ Sex) +
  geom_line(aes(colour = Year), linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(colour = Year, shape = Dataset), size = 11)+ 
  scale_color_manual(values = c("#B72779", "#2779B7"))+
  scale_shape_manual(values = c(19, 10, 8, 2, 11, 7, 46)) +
  theme_graphs()+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs()
#labs(title = "Age-Specific Mortality Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates") 
ggsave(file="Graphs/Socsim_Exp2_ASMR.jpeg", width=17, height=15, dpi=200)


#----------------------------------------------------------------------------------------------------
## Final plot combining ASFR and ASMR ----

# Years to plot limited to  two years
yrs_plot <- c("[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_whole2$age)

## Plotting ASFR and ASMR (for females) from whole SOCSIM simulation and subsets of direct ancestors and all collateral kin
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR),
          asfr_dir_wod2 %>% rename(Estimate = ASFR),
          asfr_anc_zaukgausc2 %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx),
            asmr_dir_wod2 %>% rename(Estimate = mx),
            asmr_anc_zaukgausc2 %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels),
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates"), 
         Dataset = factor(Dataset, levels =  c("Direct Ancestors (DA)", "Direct Ancestors + Collateral Kin", "Whole Simulation"))) %>%
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(aes(colour = Year), linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(colour = Year, shape = Dataset), size = 11)+ 
  scale_color_manual(values = c("#B72779", "#2779B7"))+
  scale_shape_manual(values = c(19, 18, 46)) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(x = "Age") +
  guides(shape = guide_legend(order = 1), col = guide_legend(order = 2))
ggsave(file="Graphs/Final_Socsim_Exp2_ASFR_ASMR.jpeg", width=17, height=9, dpi=200)


## For appendix
## Plotting ASFR and ASMR (for females) from whole SOCSIM simulation and subsets of direct ancestors and different collateral kin
yrs_plot_1 <- c("[1900,1905)") 

bind_rows(asfr_whole2 %>% rename(Estimate = ASFR),
          asfr_dir_wod2 %>% rename(Estimate = ASFR),
          asfr_anc_z2 %>% rename(Estimate = ASFR),
          asfr_anc_z2 %>% rename(Estimate = ASFR),
          asfr_anc_au2 %>% rename(Estimate = ASFR),
          asfr_anc_k2 %>% rename(Estimate = ASFR),
          asfr_anc_gau2 %>% rename(Estimate = ASFR),
          asfr_anc_sc2 %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx),
            asmr_dir_wod2 %>% rename(Estimate = mx),
            asmr_anc_z2 %>% rename(Estimate = mx),
            asmr_anc_z2 %>% rename(Estimate = mx),
            asmr_anc_au2 %>% rename(Estimate = mx),
            asmr_anc_k2 %>% rename(Estimate = mx),
            asmr_anc_gau2 %>% rename(Estimate = mx),
            asmr_anc_sc2 %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Dataset != "Direct Ancestors (DA)" & Year %in% yrs_plot_1) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels),
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates"), 
         Dataset = ifelse(Dataset == "DA + Siblings", "Direct Ancestors (DA) + Siblings", Dataset), 
         Dataset = factor(Dataset, levels =  c("Direct Ancestors (DA) + Siblings", "DA + Aunts/Uncles", "DA + Cousins", "DA + Great-Aunts/Uncles", 
                                               "DA + Spouse and Children", "Whole Simulation"))) %>%
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset)))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(linewidth = 1.2, color = "#B72779", show.legend = T)+ 
  geom_point(aes(shape = Dataset), color = "#B72779", size = 11)+ 
  scale_shape_manual(values = c(10, 8, 2, 11, 7, 46)) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR = scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(x = "Age")
ggsave(file="Graphs/App_Socsim_Exp2_ASFR_ASMR.jpeg", width=17, height=9, dpi=200)


## Save as .svg file for poster
# ggsave(file="Graphs/Socsim_Exp2_ASFR_ASMR.svg", device = "svg", units = "in", width=15, height=8, dpi=200) 

#----------------------------------------------------------------------------------------------------
## Summary measures: TFR and e0 ----
# Here, we use the rates by 1 year age group and 1 calendar year

# Total Fertility Rate ----
# Estimate Total Fertility Rate from asfr 1x1

# Estimate age-specific fertility rates 1x1 for the subset of direct ancestors and siblings
asfr_anc_z_1 <- map_dfr(anc_z, ~ estimate_fertility_rates(opop = .x,
                                                          final_sim_year = 2022, #[Jan-Dec]
                                                          year_min = 1750, # Closed [
                                                          year_max = 2023, # Open )
                                                          year_group = 1, 
                                                          age_min_fert = 10, # Closed [
                                                          age_max_fert = 55, # Open )
                                                          age_group = 1), # [,)
                        .id = "Sim_id") 
save(asfr_anc_z_1, file = "Measures/asfr_anc_z_1.RData")

# Estimate age-specific fertility rates 1x1 for the subset of direct ancestors and aunts/uncles
asfr_anc_au_1 <- map_dfr(anc_au, ~ estimate_fertility_rates(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed [
                                                            year_max = 2023, # Open )
                                                            year_group = 1, 
                                                            age_min_fert = 10, # Closed [
                                                            age_max_fert = 55, # Open )
                                                            age_group = 1), # [,)
                         .id = "Sim_id") 
save(asfr_anc_au_1, file = "Measures/asfr_anc_au_1.RData")

# Estimate age-specific fertility rates 1x1 for the subset of direct ancestors and cousins
asfr_anc_k_1 <- map_dfr(anc_k, ~ estimate_fertility_rates(opop = .x,
                                                          final_sim_year = 2022, #[Jan-Dec]
                                                          year_min = 1750, # Closed [
                                                          year_max = 2023, # Open )
                                                          year_group = 1, 
                                                          age_min_fert = 10, # Closed [
                                                          age_max_fert = 55, # Open )
                                                          age_group = 1), # [,)
                        .id = "Sim_id") 
save(asfr_anc_k_1, file = "Measures/asfr_anc_k_1.RData")

# Estimate age-specific fertility rates 1x1 for the subset of direct ancestors and great-aunts/uncles
asfr_anc_gau_1 <- map_dfr(anc_gau, ~ estimate_fertility_rates(opop = .x,
                                                              final_sim_year = 2022, #[Jan-Dec]
                                                              year_min = 1750, # Closed [
                                                              year_max = 2023, # Open )
                                                              year_group = 1, 
                                                              age_min_fert = 10, # Closed [
                                                              age_max_fert = 55, # Open )
                                                              age_group = 1), # [,)
                          .id = "Sim_id") 
save(asfr_anc_gau_1, file = "Measures/asfr_anc_gau_1.RData")

# Estimate age-specific fertility rates 1x1 for the subset of direct ancestors and spouse and children
asfr_anc_sc_1 <- map_dfr(anc_sc, ~ estimate_fertility_rates(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed [
                                                            year_max = 2023, # Open )
                                                            year_group = 1, 
                                                            age_min_fert = 10, # Closed [
                                                            age_max_fert = 55, # Open )
                                                            age_group = 1), # [,)
                         .id = "Sim_id") 
save(asfr_anc_sc_1, file = "Measures/asfr_anc_sc_1.RData")

# Estimate age-specific fertility rates 1x1 for the subset of direct ancestors, siblings, aunts/uncles, cousins, great-aunts/uncles, 
# spouse and children
asfr_anc_zaukgausc_1 <- map_dfr(anc_zaukgausc, ~ estimate_fertility_rates(opop = .x,
                                                                          final_sim_year = 2022, #[Jan-Dec]
                                                                          year_min = 1750, # Closed [
                                                                          year_max = 2023, # Open )
                                                                          year_group = 1, 
                                                                          age_min_fert = 10, # Closed [
                                                                          age_max_fert = 55, # Open )
                                                                          age_group = 1), # [,)
                                .id = "Sim_id") 
save(asfr_anc_zaukgausc_1, file = "Measures/asfr_anc_zaukgausc_1.RData")

# # Load mean age-specific fertility rates 1x1 for the 10 simulations
load("Measures/asfr_whole_1.RData")
## Load age-specific fertility rates 1x1 of each genealogical subset for the 10 simulations
# Load asfr 1x1 for the subset of direct ancestors
load("Measures/asfr_dir_wod_1.RData")
# Load asfr 1x1 for the subset of direct ancestors and siblings
load("Measures/asfr_anc_z_1.RData")
# Load asfr 1x1 for the subset of direct ancestors and aunts/uncles
load("Measures/asfr_anc_au_1.RData")
# Load asfr 1x1 for the subset of direct ancestors and cousins
load("Measures/asfr_anc_k_1.RData")
# Load asfr 1x1 for the subset of direct ancestors and great-aunts/uncles
load("Measures/asfr_anc_gau_1.RData")
# Load asfr 1x1 for the subset of direct ancestors and spouse and children
load("Measures/asfr_anc_sc_1.RData")
# Load asfr 1x1 for the subset of direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles, spouse and children
load("Measures/asfr_anc_zaukgausc_1.RData")

# Age breaks of fertility rates. Extract all the unique numbers from the intervals 
age_breaks_fert <- unique(as.numeric(str_extract_all(asfr_whole_1$age, "\\d+", simplify = T)))

# Retrieve age_group size
age_group_fert <- unique(diff(age_breaks_fert))

# Whole SOCSIM simulation
TFR_whole <- asfr_whole_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "Whole Simulation",
         Rate = "TFR", 
         sex = "female")

# Direct ancestors without duplicates
load("Measures/asfr_dir_wod_1.RData")
TFR_dir_wod <- asfr_dir_wod_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "Direct Ancestors (DA)",
         Rate = "TFR", 
         sex = "female")

# Direct ancestors and siblings
TFR_anc_z <- asfr_anc_z_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "DA + Siblings",
         Rate = "TFR",           
         sex = "female") 

# Direct ancestors and aunts/uncles
TFR_anc_au <- asfr_anc_au_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "DA + Aunts/Uncles",
         Rate = "TFR",           
         sex = "female") 

# Direct ancestors and cousins
TFR_anc_k <- asfr_anc_k_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "DA + Cousins",
         Rate = "TFR",
         sex = "female") 

# Direct ancestors and great-aunts/uncles
TFR_anc_gau <- asfr_anc_gau_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "DA + Great-Aunts/Uncles",
         Rate = "TFR",           
         sex = "female") 

# Direct ancestors and spouse and children
TFR_anc_sc <- asfr_anc_sc_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "DA + Spouse and Children",
         Rate = "TFR",           
         sex = "female") 

# Direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles, spouse and children
TFR_anc_zaukgausc <- asfr_anc_zaukgausc_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "TFR",           
         sex = "female")

## Plot TFR from whole SOCSIM simulation and subsets of "direct" and "extended"  family trees without duplicates

bind_rows(TFR_whole, TFR_dir_wod, TFR_anc_z, TFR_anc_z, TFR_anc_au, 
          TFR_anc_k, TFR_anc_gau, TFR_anc_sc) %>% 
  filter(Year >= 1850) %>%
  mutate(Dataset = factor(Dataset, levels =  c("Direct Ancestors (DA)", "DA + Siblings", "DA + Aunts/Uncles", "DA + Cousins", "DA + Great-Aunts/Uncles", 
                                               "DA + Spouse and Children", "Whole Simulation"))) %>% 
  ggplot(aes(x = Year, y = TFR)) +
  geom_line(aes(colour = Dataset), linewidth = 1.3)+ 
  scale_color_viridis(option = "G", discrete = T, direction = -1)+
  theme_graphs() 
# labs(title = "Total Fertility Rate in Sweden (1751-2022), 
#retrieved from a SOCSIM simulation and subsets of "direct" and "extended"  family trees without duplicates") 
ggsave(file="Graphs/socsim_Exp2_TFR.jpeg", width=17, height=9, dpi=200)

# Life Expectancy at birth ----
# Estimate life expectancy at birth from asmr 1x1 for the different genealogical subsets

# Estimate age-specific mortality rates 1x1 for the subset of direct ancestors and siblings
asmr_anc_z_1 <- map_dfr(anc_z, ~ estimate_mortality_rates(opop = .x,
                                                          final_sim_year = 2022, #[Jan-Dec]
                                                          year_min = 1750, # Closed
                                                          year_max = 2023, # Open )
                                                          year_group = 1,
                                                          age_max_mort = 110, # Open )
                                                          age_group = 1), # [,)
                        .id = "Sim_id") 
save(asmr_anc_z_1, file = "Measures/asmr_anc_z_1.RData")

# Calculate the mean and compute the life table for the subset of direct ancestors and siblings
asmr_anc_z_1 <- asmr_anc_z_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
lt_anc_z <- lt_socsim(asmr_anc_z_1)
save(lt_anc_z, file = "Measures/lt_anc_z.RData")

# Estimate age-specific mortality rates 1x1 for the subset of direct ancestors and aunts/uncles
asmr_anc_au_1 <- map_dfr(anc_au, ~ estimate_mortality_rates(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed
                                                            year_max = 2023, # Open )
                                                            year_group = 1,
                                                            age_max_mort = 110, # Open )
                                                            age_group = 1), # [,)
                         .id = "Sim_id") 
save(asmr_anc_au_1, file = "Measures/asmr_anc_au_1.RData")

# Calculate the mean and compute the life table for the subset of direct ancestors and aunts/uncles
asmr_anc_au_1 <- asmr_anc_au_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
lt_anc_au <- lt_socsim(asmr_anc_au_1)
save(lt_anc_au, file = "Measures/lt_anc_au.RData")

# Estimate age-specific mortality rates 1x1 for the subset of direct ancestors and cousins
asmr_anc_k_1 <- map_dfr(anc_k, ~ estimate_mortality_rates(opop = .x,
                                                          final_sim_year = 2022, #[Jan-Dec]
                                                          year_min = 1750, # Closed
                                                          year_max = 2023, # Open )
                                                          year_group = 1,
                                                          age_max_mort = 110, # Open )
                                                          age_group = 1), # [,)
                        .id = "Sim_id") 
save(asmr_anc_k_1, file = "Measures/asmr_anc_k_1.RData")

# Calculate the mean and compute the life table for the subset of direct ancestors and cousins 
asmr_anc_k_1 <- asmr_anc_k_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
lt_anc_k <- lt_socsim(asmr_anc_k_1)
save(lt_anc_k, file = "Measures/lt_anc_k.RData")

# Estimate age-specific mortality rates 1x1 for the subset of direct ancestors and great-aunts/uncles
asmr_anc_gau_1 <- map_dfr(anc_gau, ~ estimate_mortality_rates(opop = .x,
                                                              final_sim_year = 2022, #[Jan-Dec]
                                                              year_min = 1750, # Closed
                                                              year_max = 2023, # Open )
                                                              year_group = 1,
                                                              age_max_mort = 110, # Open )
                                                              age_group = 1), # [,)
                          .id = "Sim_id") 
save(asmr_anc_gau_1, file = "Measures/asmr_anc_gau_1.RData")

# Calculate the mean and compute the life table for the subset of direct ancestors and great-aunts/uncles, 
asmr_anc_gau_1 <- asmr_anc_gau_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
lt_anc_gau <- lt_socsim(asmr_anc_gau_1)
save(lt_anc_gau, file = "Measures/lt_anc_gau.RData")

# Estimate age-specific mortality rates 1x1 for the subset of direct ancestors, spouse and children
asmr_anc_sc_1 <- map_dfr(anc_sc, ~ estimate_mortality_rates(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed
                                                            year_max = 2023, # Open )
                                                            year_group = 1,
                                                            age_max_mort = 110, # Open )
                                                            age_group = 1), # [,)
                         .id = "Sim_id") 
save(asmr_anc_sc_1, file = "Measures/asmr_anc_sc_1.RData")

# Calculate the mean and compute the life table for the subset of direct ancestors, spouse and children
asmr_anc_sc_1 <- asmr_anc_sc_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
lt_anc_sc <- lt_socsim(asmr_anc_sc_1)
save(lt_anc_sc, file = "Measures/lt_anc_sc.RData")


# Estimate age-specific mortality rates 1x1 for the subset of direct ancestors, siblings, aunts/uncles, cousins, great-aunts/uncles, 
# spouse and children
asmr_anc_zaukgausc_1 <- map_dfr(anc_zaukgausc, ~ estimate_mortality_rates(opop = .x,
                                                                          final_sim_year = 2022, #[Jan-Dec]
                                                                          year_min = 1750, # Closed
                                                                          year_max = 2023, # Open )
                                                                          year_group = 1,
                                                                          age_max_mort = 110, # Open )
                                                                          age_group = 1), # [,)
                                .id = "Sim_id") 
save(asmr_anc_zaukgausc_1, file = "Measures/asmr_anc_zaukgausc_1.RData")

# Calculate the mean and compute the life table for the subset of direct ancestors, siblings, aunts/uncles, cousins, great-aunts/uncles, 
# spouse and children
asmr_anc_zaukgausc_1 <- asmr_anc_zaukgausc_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
lt_anc_zaukgausc <- lt_socsim(asmr_anc_zaukgausc_1)
save(lt_anc_zaukgausc, file = "Measures/lt_anc_zaukgausc.RData")


# Wrangle data for plotting

# Load life tables for the mean asmr from subset of direct ancestors
load("Measures/lt_dir_wod.RData")
# Load life tables for the mean asmr from subset of direct ancestors and siblings
load("Measures/lt_anc_z.RData")
# Load life tables for the mean asmr from subset of direct ancestors and aunts/uncles
load("Measures/lt_anc_au.RData")
# Load life tables for the mean asmr from subset of direct ancestors and cousins
load("Measures/lt_anc_k.RData")
# Load life tables for the mean asmr from subset of direct ancestors and great-aunts/uncles
load("Measures/lt_anc_gau.RData")
# Load life tables for the mean asmr from subset of direct ancestors and spouse and children
load("Measures/lt_anc_sc.RData")
# Load life tables for the mean asmr from subset of direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles, spouse and children
load("Measures/lt_anc_zaukgausc.RData")

# # Load mean age-specific mortality rates 1x1 for the 10 simulations
load("Measures/asmr_whole_1.RData")

# Year breaks. Extract all the unique numbers from the intervals 
year_breaks_mort_1 <- unique(as.numeric(str_extract_all(asmr_whole_1$year, "\\d+", simplify = T)))

# Year range to filter HMD data
year_range_mort_1 <- min(year_breaks_mort_1):max(year_breaks_mort_1-1)

# Whole SOCSIM simulation
load("Measures/lt_whole.RData")
lt_whole2 <- lt_whole %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Whole Simulation",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Direct ancestors without duplicates
load("Measures/lt_dir_wod.RData")
lt_dir_wod2 <- lt_dir_wod %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct Ancestors (DA)",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Direct ancestors and siblings
lt_anc_z2 <- lt_anc_z %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "DA + Siblings",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Direct ancestors and aunts/uncles
lt_anc_au2 <- lt_anc_au %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "DA + Aunts/Uncles",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Direct ancestors and cousins
lt_anc_k2 <- lt_anc_k %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "DA + Cousins",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Direct ancestors and great-aunts/uncles
lt_anc_gau2 <- lt_anc_gau %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "DA + Great-Aunts/Uncles",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Direct ancestors and spouse and children
lt_anc_sc2 <- lt_anc_sc %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "DA + Spouse and Children",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Direct ancestors, siblings, aunts/uncles, first cousins, great-aunt/uncles, spouse and children
lt_anc_zaukgausc2 <- lt_anc_zaukgausc %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)


bind_rows(lt_whole2, lt_dir_wod2, lt_anc_z2, lt_anc_z2, lt_anc_au2, 
          lt_anc_k2, lt_anc_gau2, lt_anc_sc2)  %>% 
  filter(Age == 0 & Year >=1850) %>%
  mutate(Sex = ifelse(sex == "female", "Female", "Male"), 
         Dataset = factor(Dataset, levels =  c("Direct Ancestors (DA)", "DA + Siblings", "DA + Aunts/Uncles", "DA + Cousins", "DA + Great-Aunts/Uncles", 
                                               "DA + Spouse and Children", "Whole Simulation"))) %>% 
  ggplot(aes(x = Year, y = ex, group = Dataset))+
  geom_line(aes(colour = Dataset), linewidth = 1.3)+
  scale_color_viridis(option = "D", discrete = T, direction = -1)+
  facet_wrap(~Sex) +
  theme_graphs() +
  labs(y = "e0") 
#  labs(title = "Life expectancy at birth in Sweden (e0), 1751-2020, retrieved from a SOCSIM simulation 
# subsets of "extended" and direct" family trees without duplicates")
ggsave(file="Graphs/socsim_Exp2_e0.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

yrs_plot2 <- c(1850, 1900, 1950, 2000)

## TFR and e0 (for females) from whole SOCSIM simulation and subsets of "direct" and all collateral kin
bind_rows(TFR_whole %>% rename(Estimate = TFR),
          TFR_dir_wod %>% rename(Estimate = TFR),
          TFR_anc_zaukgausc %>% rename(Estimate = TFR)) %>%  
  mutate(sex = "female") %>%  
  bind_rows(lt_whole2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_dir_wod2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_zaukgausc2 %>% rename(Estimate = ex) %>% filter(Age == 0)) %>%
  filter(sex == "female" & Year >= 1850) %>%
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life Expectancy at Birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life Expectancy at Birth")), 
         Dataset = factor(Dataset, levels =  c("Direct Ancestors (DA)", "Direct Ancestors + Collateral Kin", "Whole Simulation"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), 
             aes(shape = Dataset), size = 11) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#771A30", "#1A3077", "#1A7761"))+
  scale_shape_manual(values = c(19, 18, 46)) +
  #scale_x_continuous(breaks = yrs_plot)+
  theme_graphs()
# labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2020), retrieved 
# from SOCSIM microsimulation and subsets of "direct" and extended" family trees") + 
# Save the plot
ggsave(file="Graphs/Final_Socsim_Exp2_TFR_e0.jpeg", width=17, height=9, dpi=200)


## For appendix
## TFR and e0 (for females) from whole SOCSIM simulation and subsets of "direct" and different collateral kin
bind_rows(TFR_whole %>% rename(Estimate = TFR),
          TFR_dir_wod %>% rename(Estimate = TFR),
          TFR_anc_z %>% rename(Estimate = TFR),
          TFR_anc_au %>% rename(Estimate = TFR),
          TFR_anc_k %>% rename(Estimate = TFR),
          TFR_anc_gau %>% rename(Estimate = TFR),
          TFR_anc_sc %>% rename(Estimate = TFR)) %>%  
  mutate(sex = "female") %>%  
  bind_rows(lt_whole2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_dir_wod2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_z2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_au2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_k2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_gau2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_sc2 %>% rename(Estimate = ex) %>% filter(Age == 0)) %>%
  filter(sex == "female" & Dataset != "Direct Ancestors (DA)" & Year >= 1850) %>%
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life Expectancy at Birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life Expectancy at Birth")), 
         Dataset = ifelse(Dataset == "DA + Siblings", "Direct Ancestors (DA) + Siblings", Dataset), 
         Dataset = factor(Dataset, levels =  c("Direct Ancestors (DA) + Siblings", "DA + Aunts/Uncles", "DA + Cousins", "DA + Great-Aunts/Uncles", 
                                               "DA + Spouse and Children", "Whole Simulation"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), 
             aes(shape = Dataset), size = 11) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis(option = "D", discrete = T, direction =-1)+
  scale_shape_manual(values = c(10, 8, 2, 11, 7, 46)) +
  #scale_x_continuous(breaks = yrs_plot)+
  theme_graphs()
# labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2020), retrieved 
# from SOCSIM microsimulation and subsets of "direct" and extended" family trees") + 
# Save the plot
ggsave(file="Graphs/App_Socsim_Exp2_TFR_e0.jpeg", width=17, height=9, dpi=200)