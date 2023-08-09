#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Sweden - Trace and compare subsets with randomly removed proportions of children
# U:/SOCSIM/SOCSIM_Genealogies/5_Compare_Omitted_Children.R

# Remove a proportion of sub-populations (early deceased children)
# from SOCSIM microsimulations for Sweden (1751-2022) 
# Trace genealogies and compare demographic measures from the whole simulation and the subsets 

# Created on 27-06-2023
# Last modified on 04-08-2023

## NB: To run this code, it is necessary to have already run the scripts 
# 1_Run_Simulations.R, 3_Compare_Ancestors.R and 4_Compare_Kin.R
#------------------------------------------------------------------------------------------------------
## General settings and functions ----

# Prevent scientific notation (useful for the rate calculation)
options(scipen=999999)

## Load packages 
library(tidyverse)
library(ggh4x)  # To facet scales-
library(patchwork) # To combine ggplots
library(rsocsim) # Functions to estimate rates
library(svglite) # To save svg files
library(viridis)

## Load theme for the graphs and to convert SOCSIM time
source("Functions/Functions_Graphs.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with asmr calculated with rsocsim::estimate_mortality_rates()
source("Functions/Functions_Life_Table.R")

#------------------------------------------------------------------------------------------------------
## Load necessary data and randomly removed proportions of early deceased children  ----

# Load the data frame with the ancestors and relatives of 10 simulations samples
load("Subsets/anc_kin_10.RData")

# NB: This data still contains duplicates within each simulation, as someone can be a relative of different egos
# As we do not distinguish now by kin type, we will keep only unique pids before merging. 
anc_kin_10 <- anc_kin_10 %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup()  
anc_kin_list <- anc_kin_10 %>% split(.$Sim_id)

# Function to get a sample of children who died before age 1 or 5, filtered after 1751
sample_children <- function(opop = opop, age_threshold, percentage) {
  
  dead_children <- opop %>% 
    filter(dod != 0 & dod > 2400) %>% 
    mutate(age_death_months = dod-dob) %>%  
    filter(age_death_months < age_threshold*12) %>% 
    pull(pid)
  sample_size <- round(length(dead_children)*percentage/100)
  children_samp <-  data.frame(pid = sample(dead_children, sample_size, replace = F))
  
  return(children_samp)
}

# Get samples of children dead before age 1 from the genealogical subset of each simulation, using different proportions of omission
# and remove them from the anc_kin opop (i.e. the genealogical subsets.

miss_children_1_25 <- map_dfr(anc_kin_list, ~ sample_children(opop = .x,
                                                              age_threshold = 1, 
                                                              percentage = 25),
                              .id = "Sim_id") 
less_children_1_25 <- anti_join(anc_kin_10, miss_children_1_25)
save(less_children_1_25, file = "Subsets/less_children_1_25.RData")

miss_children_1_50 <- map_dfr(anc_kin_list, ~ sample_children(opop = .x,
                                                              age_threshold = 1, 
                                                              percentage = 50),
                              .id = "Sim_id")
less_children_1_50 <- anti_join(anc_kin_10, miss_children_1_50)
save(less_children_1_50, file = "Subsets/less_children_1_50.RData")

miss_children_1_75 <- map_dfr(anc_kin_list, ~ sample_children(opop = .x,
                                                              age_threshold = 5, 
                                                              percentage = 75),
                              .id = "Sim_id") 
less_children_1_75 <- anti_join(anc_kin_10, miss_children_1_75) 
save(less_children_1_75, file = "Subsets/less_children_1_75.RData")

miss_children_1_100 <- map_dfr(anc_kin_list, ~ sample_children(opop = .x,
                                                              age_threshold = 5, 
                                                              percentage = 100),
                              .id = "Sim_id") 
less_children_1_100 <- anti_join(anc_kin_10, miss_children_1_100) 
save(less_children_1_100, file = "Subsets/less_children_1_100.RData")

# Get samples of children dead before age 5 for each simulation, using different proportions of omission
miss_children_5_25 <- map_dfr(anc_kin_list, ~ sample_children(opop = .x,
                                                              age_threshold = 5, 
                                                              percentage = 25),
                              .id = "Sim_id") 
less_children_5_25 <- anti_join(anc_kin_10, miss_children_5_25) 
save(less_children_5_25, file = "Subsets/less_children_5_25.RData")

miss_children_5_50 <- map_dfr(anc_kin_list, ~ sample_children(opop = .x,
                                                              age_threshold = 5, 
                                                              percentage = 50),
                              .id = "Sim_id") 
less_children_5_50 <- anti_join(anc_kin_10, miss_children_5_50)
save(less_children_5_50, file = "Subsets/less_children_5_50.RData")


miss_children_5_75 <- map_dfr(anc_kin_list, ~ sample_children(opop = .x,
                                                              age_threshold = 5, 
                                                              percentage = 75),
                              .id = "Sim_id") 
less_children_5_75 <- anti_join(anc_kin_10, miss_children_5_75)
save(less_children_5_75, file = "Subsets/less_children_5_75.RData")

miss_children_5_100 <- map_dfr(anc_kin_list, ~ sample_children(opop = .x,
                                                              age_threshold = 5, 
                                                              percentage = 100),
                              .id = "Sim_id") 
less_children_5_100 <- anti_join(anc_kin_10, miss_children_5_100)
save(less_children_5_100, file = "Subsets/less_children_5_100.RData")

#----------------------------------------------------------------------------------------------------
## Age-Specific Fertility and Mortality rates, 5x5  -----
#  Estimate ASFR and ASMR for the genealogical subsets with direct ancestors and collateral kin, 
# after removing a proportion of children who died before age 1 or 5

# All direct ancestors and kin from the subset without 25% children dead below age 1
# Create a list of data frames with opop of the filter data
less_children_1_25 <- less_children_1_25 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 25% children dead below age 1
asfr_less_children_1_25 <- map_dfr(less_children_1_25, ~ estimate_fertility_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_min_fert = 10, 
                                                                                  age_max_fert = 55, 
                                                                                  age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_children_1_25, file = "Measures/asfr_less_children_1_25.RData")

# Estimate age-specific mortality rates from the subset without 25% children dead below age 1
asmr_less_children_1_25 <- map_dfr(less_children_1_25, ~ estimate_mortality_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_max_mort = 110, 
                                                                                  age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_children_1_25, file = "Measures/asmr_less_children_1_25.RData")


# All direct ancestors and kin from the subset without 50% children dead below age 1
# Create a list of data frames with opop of the filter data
less_children_1_50 <- less_children_1_50 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 50% children dead below age 1
asfr_less_children_1_50 <- map_dfr(less_children_1_50, ~ estimate_fertility_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_min_fert = 10, 
                                                                                  age_max_fert = 55, 
                                                                                  age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_children_1_50, file = "Measures/asfr_less_children_1_50.RData")

# Estimate age-specific mortality rates from the subset without 50% children dead below age 1
asmr_less_children_1_50 <- map_dfr(less_children_1_50, ~ estimate_mortality_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_max_mort = 110, 
                                                                                  age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_children_1_50, file = "Measures/asmr_less_children_1_50.RData")


# All direct ancestors and kin from the subset without 75% children dead below age 1
# Create a list of data frames with opop of the filter data
less_children_1_75 <- less_children_1_75 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 75% children dead below age 1
asfr_less_children_1_75 <- map_dfr(less_children_1_75, ~ estimate_fertility_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_min_fert = 10, 
                                                                                  age_max_fert = 55, 
                                                                                  age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_children_1_75, file = "Measures/asfr_less_children_1_75.RData")

# Estimate age-specific mortality rates from the subset without 75% children dead below age 1
asmr_less_children_1_75 <- map_dfr(less_children_1_75, ~ estimate_mortality_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_max_mort = 110, 
                                                                                  age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_children_1_75, file = "Measures/asmr_less_children_1_75.RData")


# All direct ancestors and kin from the subset without 100% children dead below age 1
# Create a list of data frames with opop of the filter data
less_children_1_100 <- less_children_1_100 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 75% children dead below age 1
asfr_less_children_1_100 <- map_dfr(less_children_1_100, ~ estimate_fertility_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_min_fert = 10, 
                                                                                  age_max_fert = 55, 
                                                                                  age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_children_1_100, file = "Measures/asfr_less_children_1_100.RData")

# Estimate age-specific mortality rates from the subset without 100% children dead below age 1
asmr_less_children_1_100 <- map_dfr(less_children_1_100, ~ estimate_mortality_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_max_mort = 110, 
                                                                                  age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_children_1_100, file = "Measures/asmr_less_children_1_100.RData")



# All direct ancestors and kin from the subset without 25% children dead below age 5
# Create a list of data frames with opop of the filter data
less_children_5_25 <- less_children_5_25 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 25% children dead below age 5
asfr_less_children_5_25 <- map_dfr(less_children_5_25, ~ estimate_fertility_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_min_fert = 10, 
                                                                                  age_max_fert = 55, 
                                                                                  age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_children_5_25, file = "Measures/asfr_less_children_5_25.RData")

# Estimate age-specific mortality rates from the subset without 25% children dead below age 5
asmr_less_children_5_25 <- map_dfr(less_children_5_25, ~ estimate_mortality_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_max_mort = 110, 
                                                                                  age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_children_5_25, file = "Measures/asmr_less_children_5_25.RData")


# All direct ancestors and kin from the subset without 50% children dead below age 5
# Create a list of data frames with opop of the filter data
less_children_5_50 <- less_children_5_50 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 50% children dead below age 5
asfr_less_children_5_50 <- map_dfr(less_children_5_50, ~ estimate_fertility_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_min_fert = 10, 
                                                                                  age_max_fert = 55, 
                                                                                  age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_children_5_50, file = "Measures/asfr_less_children_5_50.RData")

# Estimate age-specific mortality rates from the subset without 50% children dead below age 5
asmr_less_children_5_50 <- map_dfr(less_children_5_50, ~ estimate_mortality_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_max_mort = 110, 
                                                                                  age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_children_5_50, file = "Measures/asmr_less_children_5_50.RData")


# All direct ancestors and kin from the subset without 75% children dead below age 5
# Create a list of data frames with opop of the filter data
less_children_5_75 <- less_children_5_75 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 75% children dead below age 5
asfr_less_children_5_75 <- map_dfr(less_children_5_75, ~ estimate_fertility_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_min_fert = 10, 
                                                                                  age_max_fert = 55, 
                                                                                  age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_children_5_75, file = "Measures/asfr_less_children_5_75.RData")

# Estimate age-specific mortality rates from the subset without 75% children dead below age 5
asmr_less_children_5_75 <- map_dfr(less_children_5_75, ~ estimate_mortality_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_max_mort = 110, 
                                                                                  age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_children_5_75, file = "Measures/asmr_less_children_5_75.RData")


# All direct ancestors and kin from the subset without 100% children dead below age 5
# Create a list of data frames with opop of the filter data
less_children_5_100 <- less_children_5_100 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 100% children dead below age 5
asfr_less_children_5_100 <- map_dfr(less_children_5_100, ~ estimate_fertility_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_min_fert = 10, 
                                                                                  age_max_fert = 55, 
                                                                                  age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_children_5_100, file = "Measures/asfr_less_children_5_100.RData")

# Estimate age-specific mortality rates from the subset without 100% children dead below age 5
asmr_less_children_5_100 <- map_dfr(less_children_5_100, ~ estimate_mortality_rates(opop = .x,
                                                                                  final_sim_year = 2022, 
                                                                                  year_min = 1750, 
                                                                                  year_max = 2020, 
                                                                                  year_group = 5,
                                                                                  age_max_mort = 110, 
                                                                                  age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_children_5_100, file = "Measures/asmr_less_children_5_100.RData")

#----------------------------------------------------------------------------------------------------
## Comparison of whole simulation with genealogical subsets with different proportions of missing children ----
# ASFR ----

# Load ASFR 5x5 from the 10 simulations
load("Measures/asfr_10.RData")
# Load asfr from the subset of all direct ancestors and collateral kin
load("Measures/asfr_anc_col.RData")
# Load asfr for the genealogical subset without 25% children dead below age 1
load("Measures/asfr_less_children_1_25.RData")
# Load asfr for the genealogical subset without 50% children dead below age 1
load("Measures/asfr_less_children_1_50.RData")
# Load asfr for the genealogical subset without 75% children dead below age 1
load("Measures/asfr_less_children_1_75.RData")
# Load asfr for the genealogical subset without 100% children dead below age 1
load("Measures/asfr_less_children_1_100.RData")
# Load asfr for the genealogical subset without 25% children dead below age 5
load("Measures/asfr_less_children_5_25.RData")
# Load asfr for the genealogical subset without 50% children dead below age 5
load("Measures/asfr_less_children_5_50.RData")
# Load asfr for the genealogical subset without 75% children dead below age 5
load("Measures/asfr_less_children_5_75.RData")
# Load asfr for the genealogical subset without 100% children dead below age 5
load("Measures/asfr_less_children_5_100.RData")

## Calculate the mean of the different simulations and add relevant columns

# Whole SOCSIM simulations
asfr_whole2 <- asfr_10 %>%
  group_by(year, age) %>% 
  summarise(ASFR = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Dataset = "Whole Simulation", 
         Rate = "ASFR", 
         Omitted = NA) 

# All Direct ancestors and collateral kin 
asfr_anc_col2 <- asfr_anc_col %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "ASFR", 
         Omitted = NA) 

# All Direct ancestors and collateral kin without 25% children dead below age 1
asfr_less_children_1_25b <- asfr_less_children_1_25 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "25% Omission",
         Rate = "ASFR", 
         Omitted = "Below 1")

# All Direct ancestors and collateral kin without 50% children dead below age 1
asfr_less_children_1_50b <- asfr_less_children_1_50 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "50% Omission",
         Rate = "ASFR", 
         Omitted = "Below 1")

# All Direct ancestors and collateral kin without 75% children dead below age 1
asfr_less_children_1_75b <- asfr_less_children_1_75 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "75% Omission",
         Rate = "ASFR",
         Omitted = "Below 1")

# All Direct ancestors and collateral kin without 100% children dead below age 1
asfr_less_children_1_100b <- asfr_less_children_1_100 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "100% Omission",
         Rate = "ASFR",
         Omitted = "Below 1")


# All Direct ancestors and collateral kin without 25% children dead below age 5
asfr_less_children_5_25b <- asfr_less_children_5_25 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "25% Omission",
         Rate = "ASFR",
         Omitted = "Below 5")

# All Direct ancestors and collateral kin without 50% children dead below age 5
asfr_less_children_5_50b <- asfr_less_children_5_50 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "50% Omission",
         Rate = "ASFR",
         Omitted = "Below 5")

# All Direct ancestors and collateral kin without 75% children dead below age 5
asfr_less_children_5_75b <- asfr_less_children_5_75 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "75% Omission",
         Rate = "ASFR",
         Omitted = "Below 5")

# All Direct ancestors and collateral kin without 75% children dead below age 5
asfr_less_children_5_100b <- asfr_less_children_5_100 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "100% Omission",
         Rate = "ASFR",
         Omitted = "Below 5")

## Plot ASFR from whole simulation with genealogical subsets with different proportions of missing children

# Same years to plot. Change if necessary
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2, asfr_less_children_1_25b, asfr_less_children_1_50b, asfr_less_children_1_75b, asfr_less_children_1_100b, 
          asfr_less_children_5_25b, asfr_less_children_5_50b, asfr_less_children_5_75b,  asfr_less_children_5_100b) %>% 
  filter(year %in% yrs_plot & !is.nan(ASFR) & !is.na(Omitted)) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset, Omitted), colour = year))+
  facet_wrap(~ Omitted)+
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(shape = Dataset), size = 6)+ 
  scale_color_manual(values = c("#79B727","#B72779", "#2779B7"))+ 
  scale_shape_manual(values = c(8, 23, 22, 21, 46, 8, 23, 22, 21, 46)) +
  theme_graphs()
# labs(title = "Age-Specific Fertility Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and subsets of "direct" and "extended" family trees") 
ggsave(file="Graphs/Socsim_Exp3A_ASFR.jpeg", width=17, height=9, dpi=200)


# ASMR ----

# Load ASMR 5x5 from the 10 simulations
load("Measures/asmr_10.RData")
# Load asmr from the subset of all direct ancestors and collateral kin
load("Measures/asmr_anc_col.RData")
# Load asmr for the genealogical subset without 25% children dead below age 1
load("Measures/asmr_less_children_1_25.RData")
# Load asmr for the genealogical subset without 50% children dead below age 1
load("Measures/asmr_less_children_1_50.RData")
# Load asmr for the genealogical subset without 75% children dead below age 1
load("Measures/asmr_less_children_1_75.RData")
# Load asmr for the genealogical subset without 100% children dead below age 1
load("Measures/asmr_less_children_1_100.RData")
# Load asmr for the genealogical subset without 25% children dead below age 5
load("Measures/asmr_less_children_5_25.RData")
# Load asmr for the genealogical subset without 50% children dead below age 5
load("Measures/asmr_less_children_5_50.RData")
# Load asmr for the genealogical subset without 75% children dead below age 5
load("Measures/asmr_less_children_5_75.RData")
# Load asmr for the genealogical subset without 100% children dead below age 5
load("Measures/asmr_less_children_5_100.RData")

# Whole SOCSIM simulations
asmr_whole2 <- asmr_10 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "Whole Simulation", 
         Rate = "ASMR", 
         Omitted = NA) %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All Direct ancestors and collateral kin 
asmr_anc_col2 <- asmr_anc_col %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "ASMR", 
         Omitted = NA) %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate, Omitted)

# All Direct ancestors and collateral kin without 25% children dead below age 1
asmr_less_children_1_25b <- asmr_less_children_1_25 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "25% Omission",
         Rate = "ASMR", 
         Omitted = "Below 1") %>%
  select(year, age, Sex, mx = socsim, Dataset, Rate, Omitted)

# All Direct ancestors and collateral kin without 50% children dead below age 1
asmr_less_children_1_50b <- asmr_less_children_1_50 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "50% Omission",
         Rate = "ASMR", 
         Omitted = "Below 1") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate, Omitted)

# All Direct ancestors and collateral kin without 75% children dead below age 1
asmr_less_children_1_75b <- asmr_less_children_1_75 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "75% Omission",
         Rate = "ASMR",
         Omitted = "Below 1") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate, Omitted)

# All Direct ancestors and collateral kin without 100% children dead below age 1
asmr_less_children_1_100b <- asmr_less_children_1_100 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "100% Omission",
         Rate = "ASMR",
         Omitted = "Below 1") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate, Omitted)

# All Direct ancestors and collateral kin without 25% children dead below age 5
asmr_less_children_5_25b <- asmr_less_children_5_25 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "25% Omission",
         Rate = "ASMR",
         Omitted = "Below 5") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate, Omitted)

# All Direct ancestors and collateral kin without 50% children dead below age 5
asmr_less_children_5_50b <- asmr_less_children_5_50 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "50% Omission",
         Rate = "ASMR",
         Omitted = "Below 5") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate, Omitted)

# All Direct ancestors and collateral kin without 75% children dead below age 5
asmr_less_children_5_75b <- asmr_less_children_5_75 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "75% Omission",
         Rate = "ASMR",
         Omitted = "Below 5") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate, Omitted)

# All Direct ancestors and collateral kin without 100% children dead below age 5
asmr_less_children_5_100b <- asmr_less_children_5_100 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "100% Omission",
         Rate = "ASMR",
         Omitted = "Below 5") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate, Omitted)


## Plotting ASMR from whole simulation with genealogical subsets with different proportions of missing children

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asmr_whole2, asmr_less_children_1_25b, asmr_less_children_1_50b, asmr_less_children_1_75b,  asmr_less_children_1_100b, 
          asmr_less_children_5_25b, asmr_less_children_5_50b, asmr_less_children_5_75b, asmr_less_children_1_100b, ) %>% 
  rename(Year = year) %>% 
  filter(Year %in% yrs_plot) %>%
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx) & !is.na(Omitted)) %>% 
  mutate(Dataset = factor(Dataset, levels = c("100% Omission", "75% Omission", "50% Omission", "25% Omission", "Whole Simulation"))) %>% 
  ggplot(aes(x = age, y = mx, group = interaction(Year, Dataset, Omitted), colour = Year))+
  facet_grid(Sex ~ Omitted) +
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(shape = Dataset), size = 11)+ 
  scale_color_manual(values = c("#79B727","#B72779", "#2779B7"))+ 
  scale_shape_manual(values = c(8, 23, 22, 21, 46, 8, 23, 22, 21, 46)) +
  theme_graphs()+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs()
#labs(title = "Age-Specific Mortality Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and genealogical subsets with removed children") 
ggsave(file="Graphs/Socsim_Exp3A_ASMR.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Final plot combining ASFR and ASMR ----
# Here, we restrict to children below age 5 as the rates are calculated by 5-year age

# Years to plot 
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_whole2$age)

## Plotting ASFR and ASMR (for females) from whole SOCSIM simulation and subsets of direct ancestors and all collateral kin

By_Age_Exp3A <- 
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_anc_col2 %>% rename(Estimate = ASFR), 
          # asfr_less_children_1_25b %>% rename(Estimate = ASFR), 
          # asfr_less_children_1_50b %>% rename(Estimate = ASFR), 
          # asfr_less_children_1_75b %>% rename(Estimate = ASFR), 
          asfr_less_children_5_25b %>% rename(Estimate = ASFR), 
          asfr_less_children_5_50b %>% rename(Estimate = ASFR), 
          asfr_less_children_5_75b %>% rename(Estimate = ASFR),
          asfr_less_children_5_100b %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx), 
            asmr_anc_col2 %>% rename(Estimate = mx), 
            # asmr_less_children_1_25b %>% rename(Estimate = mx), 
            # asmr_less_children_1_50b %>% rename(Estimate = mx), 
            # asmr_less_children_1_75b %>% rename(Estimate = mx), 
            asmr_less_children_5_25b %>% rename(Estimate = mx), 
            asmr_less_children_5_50b %>% rename(Estimate = mx), 
            asmr_less_children_5_75b %>% rename(Estimate = mx),
            asmr_less_children_5_100b %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels),
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates"),
         Dataset = factor(Dataset, levels = c("100% Omission", "75% Omission", "50% Omission", "25% Omission",
                                              "Direct Ancestors + Collateral Kin", "Whole Simulation"))) %>%
    ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(shape = Dataset), size = 11)+ 
  scale_color_manual(values = c("#79B727","#B72779", "#2779B7"))+ 
  scale_shape_manual(values = c(8, 23, 22, 21, 18, 46)) + 
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() + 
  labs(x = "Age") +
  guides(shape = guide_legend(order = 1), col = guide_legend(order = 2)) +
  theme(legend.justification = "left", 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
By_Age_Exp3A
ggsave(file="Graphs/Final_Socsim_Exp3A_ASFR_ASMR.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Summary measures: TFR and e0 ----
# Here, we use the rates by 1 year age group and 1 calendar year

# Total Fertility Rate ----
# Estimate Total Fertility Rate from asfr 1x1 ----

# Estimate age-specific fertility rates 1x1 from the subset without 25% children dead below age 1
asfr_less_children_1_25_1 <- map_dfr(less_children_1_25, ~ estimate_fertility_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_min_fert = 10, 
                                                                                    age_max_fert = 55, 
                                                                                    age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_children_1_25_1, file = "Measures/asfr_less_children_1_25_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 50% children dead below age 1
asfr_less_children_1_50_1 <- map_dfr(less_children_1_50, ~ estimate_fertility_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_min_fert = 10, 
                                                                                    age_max_fert = 55, 
                                                                                    age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_children_1_50_1, file = "Measures/asfr_less_children_1_50_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 75% children dead below age 1
asfr_less_children_1_75_1 <- map_dfr(less_children_1_75, ~ estimate_fertility_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_min_fert = 10, 
                                                                                    age_max_fert = 55, 
                                                                                    age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_children_1_75_1, file = "Measures/asfr_less_children_1_75_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 100% children dead below age 1
asfr_less_children_1_100_1 <- map_dfr(less_children_1_100, ~ estimate_fertility_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_min_fert = 10, 
                                                                                    age_max_fert = 55, 
                                                                                    age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_children_1_100_1, file = "Measures/asfr_less_children_1_100_1.RData")


# Estimate age-specific fertility rates 1x1 from the subset without 25% children dead below age 5
asfr_less_children_5_25_1 <- map_dfr(less_children_5_25, ~ estimate_fertility_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_min_fert = 10, 
                                                                                    age_max_fert = 55, 
                                                                                    age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_children_5_25_1, file = "Measures/asfr_less_children_5_25_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 50% children dead below age 5
asfr_less_children_5_50_1 <- map_dfr(less_children_5_50, ~ estimate_fertility_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_min_fert = 10, 
                                                                                    age_max_fert = 55, 
                                                                                    age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_children_5_50_1, file = "Measures/asfr_less_children_5_50_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 75% children dead below age 5
asfr_less_children_5_75_1 <- map_dfr(less_children_5_75, ~ estimate_fertility_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_min_fert = 10, 
                                                                                    age_max_fert = 55, 
                                                                                    age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_children_5_75_1, file = "Measures/asfr_less_children_5_75_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 100% children dead below age 5
asfr_less_children_5_100_1 <- map_dfr(less_children_5_100, ~ estimate_fertility_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_min_fert = 10, 
                                                                                    age_max_fert = 55, 
                                                                                    age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_children_5_100_1, file = "Measures/asfr_less_children_5_100_1.RData")

# Load ASFR 1x1 and calculate TFR for plotting ----

# Load asfr 1x1 from the 10 simulations
load("Measures/asfr_10_1.RData")
# Load asfr 1x1 from the subset with all direct ancestors and collateral kin
load("Measures/asfr_anc_col_1.RData")
# Load asfr 1x1 from the subset without 25% children dead below age 1
load("Measures/asfr_less_children_1_25_1.RData")
# Load asfr 1x1 from the subset without 50% children dead below age 1
load("Measures/asfr_less_children_1_50_1.RData")
# Load asfr 1x1 from the subset without 75% children dead below age 1
load("Measures/asfr_less_children_1_75_1.RData")
# Load asfr 1x1 from the subset without 100% children dead below age 1
load("Measures/asfr_less_children_1_100_1.RData")
# Load asfr 1x1 from the subset without 25% children dead below age 5
load("Measures/asfr_less_children_5_25_1.RData")
# Load asfr 1x1 from the subset without 50% children dead below age 5
load("Measures/asfr_less_children_5_50_1.RData")
# Load asfr 1x1 from the subset without 75% children dead below age 5
load("Measures/asfr_less_children_5_75_1.RData")
# Load asfr 1x1 from the subset without 100% children dead below age 5
load("Measures/asfr_less_children_5_100_1.RData")

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

# All Direct ancestors and collateral kin
TFR_anc_col <- asfr_anc_col_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "TFR",           
         sex = "female") 

# All Direct ancestors and collateral kin without 25% children dead below age 1
TFR_less_children_1_25 <- asfr_less_children_1_25_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "25% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 50% children dead below age 1
TFR_less_children_1_50 <- asfr_less_children_1_50_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "50% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 75% children dead below age 1
TFR_less_children_1_75 <- asfr_less_children_1_75_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "75% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 100% children dead below age 1
TFR_less_children_1_100 <- asfr_less_children_1_100_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "100% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 25% children dead below age 5
TFR_less_children_5_25 <- asfr_less_children_5_25_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "25% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 50% children dead below age 5
TFR_less_children_5_50 <- asfr_less_children_5_50_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "50% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 75% children dead below age 5
TFR_less_children_5_75 <- asfr_less_children_5_75_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "75% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 100% children dead below age 5
TFR_less_children_5_100 <- asfr_less_children_5_100_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "100% Omission",
         Rate = "TFR",           
         sex = "female")

## Plot TFR from whole SOCSIM simulation and subsets with different proportions of omitted children

bind_rows(TFR_whole, TFR_anc_col,
          #TFR_less_children_1_25, TFR_less_children_1_50, TFR_less_children_1_75, TFR_less_children_1_100, 
          TFR_less_children_5_25, TFR_less_children_5_50, TFR_less_children_5_75, TFR_less_children_5_100) %>% 
  mutate(Dataset = factor(Dataset, levels = c("100% Omission", "75% Omission", "50% Omission", "25% Omission", 
                                              "Direct Ancestors + Collateral Kin", "Whole Simulation"))) %>% 
  group_by(Year, Dataset) %>% 
  summarise(TFR = mean(TFR, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year, y = TFR, group = Dataset, colour = Dataset)) +
  geom_line(linewidth = 1.3)+ 
  scale_color_viridis(option = "D", discrete = T, direction = -1)+
  theme_graphs() 
# labs(title = "Total Fertility Rate in Sweden (1751-2022), 
#retrieved from a SOCSIM simulation and subsets with different proportions of omitted children") 
ggsave(file="Graphs/Socsim_Exp3A_TFR.jpeg", width=17, height=9, dpi=200)


# Summary measure of error in TFR ----

# Difference in means
DiM_TFR_Exp3A <- bind_rows(TFR_whole, TFR_anc_col,
                    #TFR_less_children_1_25, TFR_less_children_1_50, TFR_less_children_1_75, TFR_less_children_1_100, 
                    TFR_less_children_5_25, TFR_less_children_5_50, TFR_less_children_5_75, TFR_less_children_5_100) %>%
  filter(Year > 1750) %>% 
  group_by(Year, Dataset) %>% 
  summarise(TFR = mean(TFR, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = Year, names_from = "Dataset", values_from = "TFR") %>% 
  pivot_longer(cols = 2:6, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation` , 
         Relative_Error = (Error/`Whole Simulation`)*100,
         Type = "DiM") %>% 
  select(-c(Genealogy,`Whole Simulation`)) 

# Mean of differences
MoD_TFR_Exp3A <- bind_rows(TFR_whole, TFR_anc_col,
                    #TFR_less_children_1_25, TFR_less_children_1_50, TFR_less_children_1_75, TFR_less_children_1_100, 
                    TFR_less_children_5_25, TFR_less_children_5_50, TFR_less_children_5_75, TFR_less_children_5_100) %>%
  filter(Year > 1750) %>% 
  pivot_wider(id_cols = c(Year, Sim_id), names_from = "Dataset", values_from = "TFR") %>% 
  pivot_longer(cols = 4:8, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation`, 
         Relative_Error = (Error/`Whole Simulation`)*100) %>% 
  group_by(Year, Dataset) %>% 
  reframe(Error = mean(Error, na.rm = T), 
          Relative_Error = mean(Relative_Error, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Type = "MoD")

# Bind both error measures and save the data frame
error_TFR_exp3A <- bind_rows(DiM_TFR_Exp3A, MoD_TFR_Exp3A)
save(error_TFR_exp3A, file = "Measures/error_TFR_exp3A.RData")

# Absolute Error  
error_TFR_exp3A %>% 
  ggplot(aes(x = Year, y = Error, group = Dataset, colour = Type)) +
  facet_wrap(. ~ Dataset)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3A_TFR_Error.jpeg", width=17, height=9, dpi=200)

# Relative Error
error_TFR_exp3A %>% 
  ggplot(aes(x = Year, y = Relative_Error, group = Dataset, colour = Type)) +
  facet_wrap(. ~ Dataset)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3A_TFR_Rel_Error.jpeg", width=17, height=9, dpi=200)

# Check minimum and maximum values of bias in TFR before 1900
error_TFR_exp3A %>% 
  filter(Year < 1900) %>% 
  #filter(Dataset == "25% Omission") %>% 
  filter(Dataset == "100% Omission") %>% 
  pull(Error) %>% 
  range()

# Check minimum and maximum values of relative bias in TFR before 1900
error_TFR_exp3A %>% 
  filter(Year < 1900) %>% 
  #filter(Dataset == "25% Omission") %>% 
  filter(Dataset == "100% Omission") %>% 
  pull(Relative_Error) %>% 
  range()

# Check mean values of bias in TFR before 1900
error_TFR_exp3A %>% 
  filter(Year < 1900) %>% 
  # filter(Dataset == "25% Omission") %>% 
  filter(Dataset == "100% Omission") %>% 
  pull(Error) %>% 
  mean()

# Check mean values of relative bias in TFR before 1900
error_TFR_exp3A %>% 
  filter(Year < 1900) %>% 
  # filter(Dataset == "25% Omission") %>% 
  filter(Dataset == "100% Omission") %>% 
  pull(Relative_Error) %>% 
  mean()

# Life Expectancy at birth ----
# Estimate life expectancy at birth from asmr 1x1 for the genealogical subsets ----
# with different proportions of omitted early deceased children

# Estimate age-specific mortality rates from the subset without 25% children dead below age 1
asmr_less_children_1_25_1 <- map_dfr(less_children_1_25, ~ estimate_mortality_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_max_mort = 110, 
                                                                                    age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_children_1_25_1, file = "Measures/asmr_less_children_1_25_1.RData")

# Compute life tables from the subset from the subset without 25% children dead below age 1
lt_less_children_1_25 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_children_1_25_1)
save(lt_less_children_1_25, file = "Measures/lt_less_children_1_25.RData")


# Estimate age-specific mortality rates from the subset without 50% children dead below age 1
asmr_less_children_1_50_1 <- map_dfr(less_children_1_50, ~ estimate_mortality_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_max_mort = 110, 
                                                                                    age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_children_1_50_1, file = "Measures/asmr_less_children_1_50_1.RData")

# Compute life tables from the subset from the subset without 50% children dead below age 1
lt_less_children_1_50 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_children_1_50_1)
save(lt_less_children_1_50, file = "Measures/lt_less_children_1_50.RData")


# Estimate age-specific mortality rates from the subset without 75% children dead below age 1
asmr_less_children_1_75_1 <- map_dfr(less_children_1_75, ~ estimate_mortality_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_max_mort = 110, 
                                                                                    age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_children_1_75_1, file = "Measures/asmr_less_children_1_75_1.RData")

# Compute life tables from the subset from the subset without 75% children dead below age 1
lt_less_children_1_75 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_children_1_75_1)
save(lt_less_children_1_75, file = "Measures/lt_less_children_1_75.RData")

# Estimate age-specific mortality rates from the subset without 100% children dead below age 1
asmr_less_children_1_100_1 <- map_dfr(less_children_1_100, ~ estimate_mortality_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_max_mort = 110, 
                                                                                    age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_children_1_100_1, file = "Measures/asmr_less_children_1_100_1.RData")

# Compute life tables from the subset from the subset without 100% children dead below age 1
lt_less_children_1_100 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_children_1_100_1)
save(lt_less_children_1_100, file = "Measures/lt_less_children_1_100.RData")


# Estimate age-specific mortality rates from the subset without 25% children dead below age 5
asmr_less_children_5_25_1 <- map_dfr(less_children_5_25, ~ estimate_mortality_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_max_mort = 110, 
                                                                                    age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_children_5_25_1, file = "Measures/asmr_less_children_5_25_1.RData")

# Compute life tables from the subset without 25% children dead below age 5
lt_less_children_5_25 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_children_5_25_1)
save(lt_less_children_5_25, file = "Measures/lt_less_children_5_25.RData")


# Estimate age-specific mortality rates from the subset without 50% children dead below age 5
asmr_less_children_5_50_1 <- map_dfr(less_children_5_50, ~ estimate_mortality_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_max_mort = 110, 
                                                                                    age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_children_5_50_1, file = "Measures/asmr_less_children_5_50_1.RData")

# Compute life tables from the subset without 50% children dead below age 5
lt_less_children_5_50 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_children_5_50_1)
save(lt_less_children_5_50, file = "Measures/lt_less_children_5_50.RData")


# Estimate age-specific mortality rates from the subset without 75% children dead below age 5
asmr_less_children_5_75_1 <- map_dfr(less_children_5_75, ~ estimate_mortality_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_max_mort = 110, 
                                                                                    age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_children_5_75_1, file = "Measures/asmr_less_children_5_75_1.RData")

# Compute life tables from the subset without 75% children dead below age 5
lt_less_children_5_75 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_children_5_75_1)
save(lt_less_children_5_75, file = "Measures/lt_less_children_5_75.RData")

# Estimate age-specific mortality rates from the subset without 100% children dead below age 5
asmr_less_children_5_100_1 <- map_dfr(less_children_5_100, ~ estimate_mortality_rates(opop = .x,
                                                                                    final_sim_year = 2022, 
                                                                                    year_min = 1750, 
                                                                                    year_max = 2023, 
                                                                                    year_group = 1,
                                                                                    age_max_mort = 110, 
                                                                                    age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_children_5_100_1, file = "Measures/asmr_less_children_5_100_1.RData")

# Compute life tables from the subset without 100% children dead below age 5
lt_less_children_5_100 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_children_5_100_1)
save(lt_less_children_5_100, file = "Measures/lt_less_children_5_100.RData")


# Load and wrangle life tables for plotting ----

# Load life tables from each whole SOCSIM simulation
load("Measures/lt_10.RData")
# Load life tables from subset of all direct ancestors and collateral kin
load("Measures/lt_anc_col.RData")
# Load life tables from subset without 25% children dead below age 1
load("Measures/lt_less_children_1_25.RData")
# Load life tables from subset without 50% children dead below age 1
load("Measures/lt_less_children_1_50.RData")
# Load life tables from subset without 75% children dead below age 1
load("Measures/lt_less_children_1_75.RData")
# Load life tables from subset without 100% children dead below age 1
load("Measures/lt_less_children_1_100.RData")
# Load life tables from subset without 25% children dead below age 5
load("Measures/lt_less_children_5_25.RData")
# Load life tables from subset without 50% children dead below age 5
load("Measures/lt_less_children_5_50.RData")
# Load life tables from subset without 75% children dead below age 5
load("Measures/lt_less_children_5_75.RData")
# Load life tables from subset without 100% children dead below age 5
load("Measures/lt_less_children_5_100.RData")

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

# All Direct ancestors and collateral kin 
lt_anc_col2 <- lt_anc_col %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "e0") %>%    
 select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 25% children dead below age 1
lt_less_children_1_25b <- lt_less_children_1_25 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "25% Omission",
         Rate = "e0") %>%
 select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 50% children dead below age 1
lt_less_children_1_50b <- lt_less_children_1_50 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "50% Omission",
         Rate = "e0") %>%    
 select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 75% children dead below age 1
lt_less_children_1_75b <- lt_less_children_1_75 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "75% Omission",
         Rate = "e0") %>%    
 select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 75% children dead below age 1
lt_less_children_1_100b <- lt_less_children_1_100 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "100% Omission",
         Rate = "e0") %>%    
 select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 25% children dead below age 5
lt_less_children_5_25b <- lt_less_children_5_25 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "25% Omission",
         Rate = "e0") %>%    
 select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 50% children dead below age 5
lt_less_children_5_50b <- lt_less_children_5_50 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "50% Omission",
         Rate = "e0") %>%    
 select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 75% children dead below age 5
lt_less_children_5_75b <- lt_less_children_5_75 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "75% Omission",
         Rate = "e0") %>%    
 select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 100% children dead below age 5
lt_less_children_5_100b <- lt_less_children_5_100 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "100% Omission",
         Rate = "e0") %>%    
 select(Year, Sim_id, ex, Dataset, Rate, sex, Age)


# Plot the estimates of life expectancy at birth
bind_rows(lt_whole2, lt_anc_col2,
          #lt_less_children_1_25b, lt_less_children_1_50b, lt_less_children_1_75b, lt_less_children_1_100b, 
          lt_less_children_5_25b, lt_less_children_5_50b, lt_less_children_5_75b, lt_less_children_5_100b) %>% 
  filter(Age == 0) %>%
  mutate(Dataset = factor(Dataset, levels = c("100% Omission", "75% Omission", "50% Omission", "25% Omission", 
                                              "Direct Ancestors + Collateral Kin", "Whole Simulation"))) %>% 
  group_by(Year, sex, Dataset) %>% 
  summarise(ex = mean(ex, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year, y = ex, group = Dataset, colour = Dataset)) +
  facet_grid(. ~ sex) +
  geom_line(linewidth = 1.3)+ 
  scale_color_viridis(option = "D", discrete = T, direction = -1)+
  theme_graphs() 
# labs(title = "Life expectancy at birth in Sweden (e0), 1751-2020, 
#retrieved from a SOCSIM simulation and subsets with different proportions of omitted children") 
ggsave(file="Graphs/Socsim_Exp3A_e0.jpeg", width=17, height=9, dpi=200)

# Summary measure of error in e0 ----

# Difference in means
DiM_e0_Exp3A <- bind_rows(lt_whole2, 
                   #lt_less_children_1_25b, lt_less_children_1_50b, lt_less_children_1_75b, lt_less_children_1_100b, 
                   lt_less_children_5_25b, lt_less_children_5_50b, lt_less_children_5_75b, lt_less_children_5_100b) %>%
  filter(Year > 1750 & Age == 0) %>% 
  group_by(Year, sex, Dataset) %>% 
  summarise(ex = mean(ex, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(Year:sex), names_from = "Dataset", values_from = "ex") %>% 
  pivot_longer(cols = 3:6, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation` , 
         Relative_Error = (Error/`Whole Simulation`)*100,
         Type = "DiM") %>% 
  select(-c(Genealogy,`Whole Simulation`)) 

# Mean of differences
MoD_e0_Exp3A <- bind_rows(lt_whole2, 
                   #lt_less_children_1_25b, lt_less_children_1_50b, lt_less_children_1_75b, lt_less_children_1_100b, 
                   lt_less_children_5_25b, lt_less_children_5_50b, lt_less_children_5_75b, lt_less_children_5_100b) %>%
  filter(Year > 1750 & Age == 0) %>% 
  pivot_wider(id_cols = c(Year, sex, Sim_id), names_from = "Dataset", values_from = "ex") %>% 
  pivot_longer(cols = 5:8, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation`, 
         Relative_Error = (Error/`Whole Simulation`)*100) %>% 
  group_by(Year, sex, Dataset) %>% 
  reframe(Error = mean(Error, na.rm = T), 
          Relative_Error = mean(Relative_Error, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Type = "MoD")

# Bind both error measures and save the data frame
error_e0_exp3A <- bind_rows(DiM_e0_Exp3A, MoD_e0_Exp3A) 
save(error_e0_exp3A, file = "Measures/error_e0_exp3A.RData")

# Absolute error
error_e0_exp3A %>%
  ggplot(aes(x = Year, y = Error, group = Dataset, colour = Type)) +
  facet_grid(Dataset ~ sex)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3A_e0_Error.jpeg", width=17, height=9, dpi=200)

# Relative error
error_e0_exp3A %>% 
  ggplot(aes(x = Year, y = Relative_Error, group = Dataset, colour = Type)) +
  facet_grid(Dataset ~ sex)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3A_e0_Rel_Error.jpeg", width=17, height=9, dpi=200)


# Check minimum and maximum values of bias in e0 over the whole period for 25% omission
error_e0_exp3A %>% 
  filter(sex == "female" & Error >0) %>%
  filter(Dataset == "25% Omission") %>% 
  #filter(Dataset == "100% Omission") %>% 
  pull(Error) %>%
  range()

# Check minimum and maximum values of bias in e0 over the whole period for 100% omission
error_e0_exp3A %>% 
  filter(sex == "female") %>%
  # filter(Dataset == "25% Omission") %>% 
  filter(Dataset == "100% Omission") %>% 
  pull(Error) %>%
  range()

# Check minimum and maximum values of relative bias in e0 over the whole period for 25% omission
error_e0_exp3A %>% 
  filter(sex == "female" & Error >0) %>%
  filter(Dataset == "25% Omission") %>% 
  #filter(Dataset == "100% Omission") %>% 
  pull(Relative_Error) %>%
  range()

# Check minimum and maximum values of relative bias in e0 over the whole period for 100% omission
error_e0_exp3A %>% 
  filter(sex == "female") %>%
  # filter(Dataset == "25% Omission") %>% 
  filter(Dataset == "100% Omission") %>% 
  pull(Relative_Error) %>%
  range()

#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

yrs_plot2 <- c(1750, 1800, 1850, 1900, 1950, 2000)

## TFR and e0 (for females) from whole SOCSIM simulation and subsets of "direct" and all collateral kin

Summary_Exp3A <- 
bind_rows(TFR_whole %>% rename(Estimate = TFR),
          TFR_anc_col %>% rename(Estimate = TFR),
          # TFR_less_children_1_25 %>% rename(Estimate = TFR), 
          # TFR_less_children_1_50 %>% rename(Estimate = TFR), 
          # TFR_less_children_1_75 %>% rename(Estimate = TFR), 
          TFR_less_children_5_25 %>% rename(Estimate = TFR), 
          TFR_less_children_5_50 %>% rename(Estimate = TFR), 
          TFR_less_children_5_75 %>% rename(Estimate = TFR),
          TFR_less_children_5_100 %>% rename(Estimate = TFR)) %>%  
  mutate(sex = "female") %>%  
  bind_rows(lt_whole2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_col2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            # lt_less_children_1_25b %>% rename(Estimate = ex) %>% filter(Age == 0), 
            # lt_less_children_1_50b %>% rename(Estimate = ex) %>% filter(Age == 0), 
            # lt_less_children_1_75b %>% rename(Estimate = ex) %>% filter(Age == 0), 
            lt_less_children_5_25b %>% rename(Estimate = ex) %>% filter(Age == 0), 
            lt_less_children_5_50b %>% rename(Estimate = ex) %>% filter(Age == 0), 
            lt_less_children_5_75b %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_less_children_5_100b %>% rename(Estimate = ex) %>% filter(Age == 0)) %>%
  filter(sex == "female") %>%
  group_by(Year, Dataset, Rate, sex) %>% 
  summarise(Estimate = mean(Estimate, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life Expectancy at Birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life Expectancy at Birth")),
         Dataset = factor(Dataset, levels = c("100% Omission", "75% Omission", "50% Omission", "25% Omission",
                                              "Direct Ancestors + Collateral Kin", "Whole Simulation"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), 
             aes(shape = Dataset), size = 11) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#FFBE4D", "#FF834C", "#E7495B", "#B90E6D","#75007A", "#007A75"))+
  scale_shape_manual(values = c(8, 23, 22, 21, 18, 46)) +
  theme_graphs() +
  theme(legend.justification = "left",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
Summary_Exp3A
# labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and subsets with different proportions of omitted children")
# Save the plot
ggsave(file="Graphs/Final_Socsim_Exp3A_TFR_e0.jpeg", width=17, height=9, dpi=200)
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

By_Age_Exp3A + 
  geom_text(data = plot_labs1, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") + 
  theme(plot.margin = margin(0,0,1,0, "cm")) +
  Summary_Exp3A + 
  geom_text(data = plot_labs2, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") +
  plot_layout(ncol = 1)
ggsave(file="Graphs/Final_Socsim_Exp3A_Combined.jpeg", width=18, height=21, dpi=200)