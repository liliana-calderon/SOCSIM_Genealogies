#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Sweden - Trace and compare subsets with randomly removed proportions of childless women
# U:/SOCSIM/SOCSIM_Genealogies/6_Compare_Omitted_women.R

# Remove a proportion of sub-populations (childless women)
# from SOCSIM microsimulations for Sweden (1751-2022) 
# Trace genealogies and compare demographic measures from the whole simulation and the subsets 

# Created on 11-07-2023
# Last modified on 03-09-2024

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
library(viridis)

## Load theme for the graphs and to convert SOCSIM time
source("Functions/Functions_Graphs.R")

## Load functions to estimate age-specific fertility and mortality rates 
# These are a slightly modified version of the functions in the rsocsim package 
# that allow to handle the intentional duplicates in the data (in fertility rates)
# and retrieve the last month from max(dod) instead of dob. 
# Important here as max(dob) in most subsets and simulations is not the last simulated month
# So, the function in rsocsim would assign incorrectly the last month of the simulation
# and hence calculate wrongly the years of birth and death
# Also, the direct ancestors offspring (e.g., ego, siblings, aunts/uncles etc) 
# are only counted in the numerator but not in the denominator of fertility rates
source("Functions/Functions_Fertility_Rates_Mod.R")
source("Functions/Functions_Mortality_Rates_Mod.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with asmr calculated with estimate_mortality_rates_mod()
source("Functions/Functions_Life_Table.R")
#------------------------------------------------------------------------------------------------------
## Load necessary data and randomly removed proportions of childless women ----

# Load the data frame with the ancestors and their offspring of 10 simulations samples
load("Subsets/anc_kin_10.RData")

# NB: This data still contains duplicates within each simulation, as someone can be a relative of different egos
# As we do not distinguish now between type of kin, we will keep only unique pids before merging. 
### COUSINS NEED TO BE REMOVED
anc_kin_10 <- anc_kin_10 %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup()  
anc_kin_list <- anc_kin_10 %>% split(.$Sim_id)


# Function to get a sample of childless women, born after 1735
# who survived at least until reproductive ages

sample_women <- function(opop = opop, percentage) {
  
  childless_women <- opop %>% 
    filter(fem == 1 & lborn == 0) %>% 
    mutate(age_death_months = dod-dob) %>%  
    # Filtered after 1735 to keep women of at least 15 years old in 1751, and who died after age 15
    filter(dob > 2220 & age_death_months > 15*12) %>% 
    pull(pid)
  sample_size <- round(length(childless_women)*percentage/100)
  women_samp <-  data.frame(pid = sample(childless_women, sample_size, replace = F))
  
  return(women_samp)
}

# Get samples of childless women from the genealogical subset of each simulation, using different proportions of omission
# and remove them from the anc_kin opop (i.e. the genealogical subsets)

miss_women_25 <- map_dfr(anc_kin_list, ~ sample_women(opop = .x, percentage = 25),
                         .id = "Sim_id") 
less_women_25 <- anti_join(anc_kin_10, miss_women_25)
save(less_women_25, file = "Subsets/less_women_25.RData")

miss_women_50 <- map_dfr(anc_kin_list, ~ sample_women(opop = .x, percentage = 50),
                         .id = "Sim_id")
less_women_50 <- anti_join(anc_kin_10, miss_women_50)
save(less_women_50, file = "Subsets/less_women_50.RData")

miss_women_75 <- map_dfr(anc_kin_list, ~ sample_women(opop = .x, percentage = 75),
                         .id = "Sim_id") 
less_women_75 <- anti_join(anc_kin_10, miss_women_75) 
save(less_women_75, file = "Subsets/less_women_75.RData")

miss_women_100 <- map_dfr(anc_kin_list, ~ sample_women(opop = .x, percentage = 100),
                          .id = "Sim_id") 
less_women_100 <- anti_join(anc_kin_10, miss_women_100) 
save(less_women_100, file = "Subsets/less_women_100.RData")

#----------------------------------------------------------------------------------------------------
## Age-Specific Fertility and Mortality rates, 5x5  -----
#  Estimate ASFR and ASMR for the genealogical subsets with direct ancestors and their offspring, 
# after removing a proportion of childless women

# All direct ancestors and kin from the subset without 25% childless women
# Create a list of data frames with opop of the filter data
less_women_25 <- less_women_25 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 25% childless women
asfr_less_women_25 <- map_dfr(less_women_25, ~ estimate_fertility_rates_mod(opop = .x,
                                                                            final_sim_year = 2072, 
                                                                            year_min = 1750, 
                                                                            year_max = 2020, 
                                                                            year_group = 5,
                                                                            age_min_fert = 10, 
                                                                            age_max_fert = 55, 
                                                                            age_group = 5), 
                              .id = "Sim_id") 
save(asfr_less_women_25, file = "Measures/asfr_less_women_25.RData")

# Estimate age-specific mortality rates from the subset without 25% childless women
asmr_less_women_25 <- map_dfr(less_women_25, ~ estimate_mortality_rates_mod(opop = .x,
                                                                            final_sim_year = 2072, 
                                                                            year_min = 1750, 
                                                                            year_max = 2020, 
                                                                            year_group = 5,
                                                                            age_max_mort = 110,
                                                                            age_group = 5),
                              .id = "Sim_id") 
save(asmr_less_women_25, file = "Measures/asmr_less_women_25.RData")


# All direct ancestors and kin from the subset without 25% childless women
# Create a list of data frames with opop of the filter data
less_women_50 <- less_women_50 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 25% childless women
asfr_less_women_50 <- map_dfr(less_women_50, ~ estimate_fertility_rates_mod(opop = .x,
                                                                            final_sim_year = 2072,
                                                                            year_min = 1750, 
                                                                            year_max = 2020, 
                                                                            year_group = 5,
                                                                            age_min_fert = 10, 
                                                                            age_max_fert = 55, 
                                                                            age_group = 5), 
                              .id = "Sim_id") 
save(asfr_less_women_50, file = "Measures/asfr_less_women_50.RData")

# Estimate age-specific mortality rates from the subset without 25% childless women
asmr_less_women_50 <- map_dfr(less_women_50, ~ estimate_mortality_rates_mod(opop = .x,
                                                                            final_sim_year = 2072,
                                                                            year_min = 1750, 
                                                                            year_max = 2020, 
                                                                            year_group = 5,
                                                                            age_max_mort = 110, 
                                                                            age_group = 5),
                              .id = "Sim_id") 
save(asmr_less_women_50, file = "Measures/asmr_less_women_50.RData")


# All direct ancestors and kin from the subset without 75% childless women
# Create a list of data frames with opop of the filter data
less_women_75 <- less_women_75 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 75% childless women
asfr_less_women_75 <- map_dfr(less_women_75, ~ estimate_fertility_rates_mod(opop = .x,
                                                                            final_sim_year = 2072, 
                                                                            year_min = 1750, 
                                                                            year_max = 2020,
                                                                            year_group = 5,
                                                                            age_min_fert = 10, 
                                                                            age_max_fert = 55, 
                                                                            age_group = 5), 
                              .id = "Sim_id") 
save(asfr_less_women_75, file = "Measures/asfr_less_women_75.RData")

# Estimate age-specific mortality rates from the subset without 75% childless women
asmr_less_women_75 <- map_dfr(less_women_75, ~ estimate_mortality_rates_mod(opop = .x,
                                                                            final_sim_year = 2072, 
                                                                            year_min = 1750, 
                                                                            year_max = 2020, 
                                                                            year_group = 5,
                                                                            age_max_mort = 110, 
                                                                            age_group = 5),
                              .id = "Sim_id") 
save(asmr_less_women_75, file = "Measures/asmr_less_women_75.RData")


# All direct ancestors and kin from the subset without 100% childless women
# Create a list of data frames with opop of the filter data
less_women_100 <- less_women_100 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 100% childless women
asfr_less_women_100 <- map_dfr(less_women_100, ~ estimate_fertility_rates_mod(opop = .x,
                                                                              final_sim_year = 2072, 
                                                                              year_min = 1750, 
                                                                              year_max = 2020,
                                                                              year_group = 5,
                                                                              age_min_fert = 10, 
                                                                              age_max_fert = 55, 
                                                                              age_group = 5), 
                               .id = "Sim_id") 
save(asfr_less_women_100, file = "Measures/asfr_less_women_100.RData")

# Estimate age-specific mortality rates from the subset without 100% childless women
asmr_less_women_100 <- map_dfr(less_women_75, ~ estimate_mortality_rates_mod(opop = .x,
                                                                             final_sim_year = 2072, 
                                                                             year_min = 1750, 
                                                                             year_max = 2020, 
                                                                             year_group = 5,
                                                                             age_max_mort = 110, 
                                                                             age_group = 5),
                               .id = "Sim_id") 
save(asmr_less_women_100, file = "Measures/asmr_less_women_100.RData")

#----------------------------------------------------------------------------------------------------
## Comparison of whole simulation with genealogical subsets with different proportions of missing childless women ----
# ASFR ----

# Load ASFR 5x5 from the 10 simulations
load("Measures/asfr_10.RData")
# Load asfr from the subset of All direct ancestors and their offspring
load("Measures/asfr_anc_off.RData")
# Load asfr for the genealogical subset without 25% childless women
load("Measures/asfr_less_women_25.RData")
# Load asfr for the genealogical subset without 25% childless women
load("Measures/asfr_less_women_50.RData")
# Load asfr for the genealogical subset without 75% childless women
load("Measures/asfr_less_women_75.RData")
# Load asfr for the genealogical subset without 100% childless women
load("Measures/asfr_less_women_100.RData")

## Calculate the mean of the different simulations and add relevant columns

# Whole SOCSIM simulations
asfr_whole2 <- asfr_10 %>%
  group_by(year, age) %>% 
  summarise(ASFR = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Dataset = "Whole Simulation", 
         Rate = "ASFR") 

# All direct ancestors and their offspring 
asfr_anc_off2 <- asfr_anc_off %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct Ancestors and their Offspring",
         Rate = "ASFR") 

# All direct ancestors and their offspring without 25% childless women
asfr_less_women_25b <- asfr_less_women_25 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "25% Omission",
         Rate = "ASFR")

# All direct ancestors and their offspring without 50% childless women
asfr_less_women_50b <- asfr_less_women_50 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "50% Omission",
         Rate = "ASFR")

# All direct ancestors and their offspring without 75% childless women
asfr_less_women_75b <- asfr_less_women_75 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "75% Omission",
         Rate = "ASFR")

# All direct ancestors and their offspring without 100% childless women
asfr_less_women_100b <- asfr_less_women_100 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "100% Omission",
         Rate = "ASFR")

## Plot ASFR from whole simulation with genealogical subsets with different proportions of missing childless women

# Same years to plot than above (in intervals). 
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2, asfr_less_women_25b, asfr_less_women_50b, asfr_less_women_75b, asfr_less_women_100b) %>% 
  filter(year %in% yrs_plot & !is.nan(ASFR)) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset), colour = year))+
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(shape = Dataset), size = 6)+ 
  scale_color_manual(values = c("#79B727", "#2779B7", "#B72779"))+ 
  scale_shape_manual(values = c(8, 21, 22, 23, 46)) +
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3B_ASFR.jpeg", width=17, height=9, dpi=300)


# ASMR ----

# Load ASMR 5x5 from the 10 simulations
load("Measures/asmr_10.RData")
# Load asmr from the subset of All direct ancestors and their offspring, calculated on 3_Compare_Ancestors
#load("Measures/asmr_anc_off.RData")
load("Measures/asmr_anc_off.RData")
# Load asmr for the genealogical subset without 25% childless women
load("Measures/asmr_less_women_25.RData")
# Load asmr for the genealogical subset without 25% childless women
load("Measures/asmr_less_women_50.RData")
# Load asmr for the genealogical subset without 75% childless women
load("Measures/asmr_less_women_75.RData")
# Load asmr for the genealogical subset without 100% childless women
load("Measures/asmr_less_women_100.RData")

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_anc_off$age)

# Whole SOCSIM simulations
asmr_whole2 <- asmr_10 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "Whole Simulation", 
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

asmr_anc_off2 <- asmr_anc_off %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "Direct Ancestors and their Offspring",
         Rate = "ASMR", 
         Omitted = "NA") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All direct ancestors and their offspring without 25% childless women
asmr_less_women_25b <- asmr_less_women_25 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "25% Omission",
         Rate = "ASMR") %>%
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All direct ancestors and their offspring without 50% childless women
asmr_less_women_50b <- asmr_less_women_50 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "50% Omission",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All direct ancestors and their offspring without 75% childless women
asmr_less_women_75b <- asmr_less_women_75 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "75% Omission",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All direct ancestors and their offspring without 75% childless women
asmr_less_women_100b <- asmr_less_women_100 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "100% Omission",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)


## Plotting ASMR from whole simulation with genealogical subsets with different proportions of missing childless women

# Same years to plot than above (in intervals). 
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asmr_whole2, asmr_less_women_25b, asmr_less_women_50b, asmr_less_women_75b, asmr_less_women_100b) %>% 
  rename(Year = year) %>% 
  filter(Year %in% yrs_plot) %>%
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>% 
  ggplot(aes(x = age, y = mx, group = interaction(Year, Dataset), colour = Year))+
  facet_grid(. ~ Sex) +
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(shape = Dataset), size = 11)+ 
  scale_color_manual(values = c("#79B727", "#2779B7", "#B72779"))+ 
  scale_shape_manual(values = c(8, 21, 22, 23, 46)) +
  theme_graphs()+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3B_ASMR.jpeg", width=17, height=9, dpi=300)

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

## Plotting ASFR and ASMR (for females) from whole SOCSIM simulation and genealogical subsets
By_Age_Exp3B <- 
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_anc_off2 %>% rename(Estimate = ASFR), 
          asfr_less_women_25b %>% rename(Estimate = ASFR),
          # asfr_less_women_50b %>% rename(Estimate = ASFR),
          # asfr_less_women_75b %>% rename(Estimate = ASFR), 
          asfr_less_women_100b %>% rename(Estimate = ASFR)) %>%
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx), 
            asmr_anc_off2 %>% rename(Estimate = mx), 
            asmr_less_women_25b %>% rename(Estimate = mx),
            # asmr_less_women_50b %>% rename(Estimate = mx),
            # asmr_less_women_75b %>% rename(Estimate = mx),
            asmr_less_women_100b %>% rename(Estimate = mx)) %>%
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot1) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels),
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates"),
         Dataset = factor(Dataset, levels = c("100% Omission", "25% Omission", 
                                              "Direct Ancestors and their Offspring", "Whole Simulation"))) %>%
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(linewidth = 1.3, show.legend = T)+
  geom_point(data = . %>% filter(age %in% age_plot), 
             aes(shape = Dataset), size = 11) +
  scale_color_manual(values = c("#2779B7"))+ 
  scale_shape_manual(values = c(8, 22, 18, 46)) + 
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() + 
  labs(x = "Age") +
  guides(colour = "none") +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
By_Age_Exp3B

## Plot ASFR and ASMR (for females), with three years for appendix

# Choose three years to plot
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels for fertility to fix them across all plots
age_levels_asfr <- levels(asfr_whole2$age)

bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_anc_off2 %>% rename(Estimate = ASFR), 
          asfr_less_women_25b %>% rename(Estimate = ASFR),
          # asfr_less_women_50b %>% rename(Estimate = ASFR),
          # asfr_less_women_75b %>% rename(Estimate = ASFR), 
          asfr_less_women_100b %>% rename(Estimate = ASFR)) %>%
  complete(year, age, Dataset, Rate, fill = list(Estimate = NA)) %>% 
  mutate(Sex = "Female") %>%
  bind_rows(asmr_whole2 %>% rename(Estimate = mx), 
            asmr_anc_off2 %>% rename(Estimate = mx), 
            asmr_less_women_25b %>% rename(Estimate = mx),
            # asmr_less_women_50b %>% rename(Estimate = mx),
            # asmr_less_women_75b %>% rename(Estimate = mx),
            asmr_less_women_100b %>% rename(Estimate = mx)) %>%
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels),
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates"),
         Dataset = factor(Dataset, levels = c("100% Omission", "25% Omission", 
                                              "Direct Ancestors and their Offspring", "Whole Simulation"))) %>%
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(Year ~ Rate, nrow = 3, ncol = 2, scales = "free") + 
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(data = . %>% filter(age %in% age_plot), 
             aes(shape = Dataset), size = 11) +
  scale_color_manual(values = c("#79B727", "#2779B7", "#B72779"))+ 
  scale_shape_manual(values = c(8, 22, 18, 46)) + 
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10"),
                               ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10"),
                               ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10")), 
                      x = list(ASFR = scale_x_discrete(limits = age_levels_asfr, guide = guide_axis(angle = 90)),
                               ASMR =  scale_x_discrete(limits = age_levels, guide = guide_axis(angle = 90)),
                               ASFR = scale_x_discrete(limits = age_levels_asfr, guide = guide_axis(angle = 90)),
                               ASMR =  scale_x_discrete(limits = age_levels, guide = guide_axis(angle = 90)),
                               ASFR = scale_x_discrete(limits = age_levels_asfr, guide = guide_axis(angle = 90)),
                               ASMR =  scale_x_discrete(limits = age_levels, guide = guide_axis(angle = 90))))+
  theme_graphs() + 
  labs(x = "Age") +
  guides(colour = "none") +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
ggsave(file="Final_Graphs/App_Socsim_Exp3B_ASFR_ASMR.jpeg", width=20, height=25, dpi=300)

#----------------------------------------------------------------------------------------------------
# Figure for EPC presentation

# Choose three years to plot
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 
age_plot <- c("[0,1)", "[1,5)", "[10,15)", "[15,20)", "[20,25)", "[30,35)", "[40,45)", "[50,55)",  "[60,65)", 
              "[70,75)", "[80,85)", "[90,95)", "[100,105)") 

bind_rows(asmr_whole2, asmr_anc_off2, asmr_less_women_25b, asmr_less_women_100b) %>%
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
  scale_shape_manual(values = c(8, 22, 18, 46)) + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs() +
  labs(x = "Age") +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  guides(colour = "none")

ggsave(file="Graphs/Socsim_Exp3B_ASMR_years.jpeg", width=24, height=9, dpi=300)
#----------------------------------------------------------------------------------------------------
## Summary measures: TFR and e0 ----
# Here, we use the rates by 1 year age group and 1 calendar year

# Total Fertility Rate ----
# Estimate Total Fertility Rate from asfr 1x1 ----

# Estimate age-specific fertility rates 1x1 from the subset without 25% childless women
asfr_less_women_25_1 <- map_dfr(less_women_25, ~ estimate_fertility_rates_mod(opop = .x,
                                                                              final_sim_year = 2072,
                                                                              year_min = 1750,
                                                                              year_max = 2023,
                                                                              year_group = 1,
                                                                              age_min_fert = 10,
                                                                              age_max_fert = 55, 
                                                                              age_group = 1), 
                                .id = "Sim_id") 
save(asfr_less_women_25_1, file = "Measures/asfr_less_women_25_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 25% childless women
asfr_less_women_50_1 <- map_dfr(less_women_50, ~ estimate_fertility_rates_mod(opop = .x,
                                                                              final_sim_year = 2072,
                                                                              year_min = 1750, 
                                                                              year_max = 2023, 
                                                                              year_group = 1,
                                                                              age_min_fert = 10, 
                                                                              age_max_fert = 55, 
                                                                              age_group = 1), 
                                .id = "Sim_id") 
save(asfr_less_women_50_1, file = "Measures/asfr_less_women_50_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 75% childless women
asfr_less_women_75_1 <- map_dfr(less_women_75, ~ estimate_fertility_rates_mod(opop = .x,
                                                                              final_sim_year = 2072, 
                                                                              year_min = 1750, 
                                                                              year_max = 2023, 
                                                                              year_group = 1,
                                                                              age_min_fert = 10, 
                                                                              age_max_fert = 55, 
                                                                              age_group = 1), 
                                .id = "Sim_id") 
save(asfr_less_women_75_1, file = "Measures/asfr_less_women_75_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 100% childless women
asfr_less_women_100_1 <- map_dfr(less_women_100, ~ estimate_fertility_rates_mod(opop = .x,
                                                                                final_sim_year = 2072, 
                                                                                year_min = 1750, 
                                                                                year_max = 2023, 
                                                                                year_group = 1,
                                                                                age_min_fert = 10, 
                                                                                age_max_fert = 55, 
                                                                                age_group = 1), 
                                 .id = "Sim_id") 
save(asfr_less_women_100_1, file = "Measures/asfr_less_women_100_1.RData")


# Load ASFR 1x1 and calculate TFR for plotting ----

# Load asfr 1x1 from the 10 simulations
load("Measures/asfr_10_1.RData")
# Load asfr 1x1 from the subset with All direct ancestors and their offspring
load("Measures/asfr_anc_off_1.RData")
# Load asfr 1x1 from the subset without 25% childless women
load("Measures/asfr_less_women_25_1.RData")
# Load asfr 1x1 from the subset without 25% childless women
load("Measures/asfr_less_women_50_1.RData")
# Load asfr 1x1 from the subset without 75% childless women
load("Measures/asfr_less_women_75_1.RData")
# Load asfr 1x1 from the subset without 100% childless women
load("Measures/asfr_less_women_100_1.RData")

# Age breaks of fertility rates. Extract all the unique numbers from the intervals 
age_breaks_fert_1 <- unique(as.numeric(str_extract_all(asfr_10_1$age, "\\d+", simplify = T)))

# Whole SOCSIM simulations
TFR_whole <- asfr_10_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim, na.rm = T)) %>%
  ungroup() %>% 
  mutate(Dataset = "Whole Simulation",
         Rate = "TFR", 
         sex = "female")

# All direct ancestors and their offspring
TFR_anc_off <- asfr_anc_off_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  # Some ages can have infinite (N_Births/0_Pop) and NaN (0_Births/0_Pop) values
  filter(!is.infinite(socsim)) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim, na.rm = T)) %>%
  ungroup() %>% 
  mutate(Dataset = "Direct Ancestors and their Offspring",
         Rate = "TFR",           
         sex = "female") 

# All direct ancestors and their offspring without 25% childless women
TFR_less_women_25 <- asfr_less_women_25_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  # Some ages can have infinite (N_Births/0_Pop) and NaN (0_Births/0_Pop) values
  filter(!is.infinite(socsim)) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim, na.rm = T)) %>%
  ungroup() %>% 
  mutate(Dataset = "25% Omission",
         Rate = "TFR",           
         sex = "female")

# All direct ancestors and their offspring without 50% childless women
TFR_less_women_50 <- asfr_less_women_50_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  # Some ages can have infinite (N_Births/0_Pop) and NaN (0_Births/0_Pop) values
  filter(!is.infinite(socsim)) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim, na.rm = T)) %>%
  ungroup() %>% 
  mutate(Dataset = "50% Omission",
         Rate = "TFR",           
         sex = "female")

# All direct ancestors and their offspring without 75% childless women
TFR_less_women_75 <- asfr_less_women_75_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  # Some ages can have infinite (N_Births/0_Pop) and NaN (0_Births/0_Pop) values
  filter(!is.infinite(socsim)) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim, na.rm = T)) %>%
  ungroup() %>% 
  mutate(Dataset = "75% Omission",
         Rate = "TFR",           
         sex = "female")

# All direct ancestors and their offspring without 100% childless women
TFR_less_women_100 <- asfr_less_women_100_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  # Some ages can have infinite (N_Births/0_Pop) and NaN (0_Births/0_Pop) values
  filter(!is.infinite(socsim)) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim, na.rm = T)) %>%
  ungroup() %>% 
  mutate(Dataset = "100% Omission",
         Rate = "TFR",           
         sex = "female")

## Plot TFR from whole SOCSIM simulation and subsets with different proportions of omitted children

bind_rows(TFR_whole, TFR_anc_off,
          TFR_less_women_25, TFR_less_women_50, TFR_less_women_75, TFR_less_women_100) %>% 
  group_by(Year, Dataset) %>% 
  summarise(TFR = mean(TFR, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year, y = TFR, group = Dataset, colour = Dataset)) +
  geom_line(linewidth = 1.3)+ 
  scale_color_viridis(option = "D", discrete = T, direction = -1)+
  theme_graphs() 
ggsave(file="Graphs/Socsim_Exp3B_TFR.jpeg", width=17, height=9, dpi=300)


# Summary measure of error in TFR ----

# Difference in means
DiM_TFR_Exp3B <- bind_rows(TFR_whole, TFR_anc_off,
           TFR_less_women_25, TFR_less_women_50, TFR_less_women_75, TFR_less_women_100) %>%
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
MoD_TFR_Exp3B <- bind_rows(TFR_whole, TFR_anc_off,
                    TFR_less_women_25, TFR_less_women_50, TFR_less_women_75, TFR_less_women_100) %>%
  filter(Year > 1750) %>% 
  pivot_wider(id_cols = c(Year, Sim_id), names_from = "Dataset", values_from = "TFR") %>% 
  pivot_longer(cols = 4:8, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation`, 
         Relative_Error = (Error/`Whole Simulation`)*100) %>% 
  group_by(Year, Dataset) %>% 
  summarise(Error = mean(Error, na.rm = T), 
          Relative_Error = mean(Relative_Error, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Type = "MoD")

# Bind both error measures and save the data frame
error_TFR_exp3B <- bind_rows(DiM_TFR_Exp3B, MoD_TFR_Exp3B)
save(error_TFR_exp3B, file = "Measures/error_TFR_exp3B.RData")

# Absolute Error  
error_TFR_exp3B %>% 
  ggplot(aes(x = Year, y = Error, group = Dataset, colour = Type)) +
  facet_wrap(. ~ Dataset)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3B_TFR_Error.jpeg", width=17, height=9, dpi=300)

# Relative Error
error_TFR_exp3B %>% 
  ggplot(aes(x = Year, y = Relative_Error, group = Dataset, colour = Type)) +
  facet_wrap(. ~ Dataset)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3B_TFR_Rel_Error.jpeg", width=17, height=9, dpi=300)

# Check minimum and maximum values of bias in TFR 
error_TFR_exp3B %>% 
  #filter(Dataset == "25% Omission") %>% # 0.1359794 1.8437241
  filter(Dataset == "100% Omission") %>% # 0.05434528 1.78261301
  pull(Error) %>% 
  range()

# Check mean values of bias in TFR
error_TFR_exp3B %>% 
  # filter(Dataset == "25% Omission") %>% # 0.4735607
  # filter(Dataset == "100% Omission") %>% # 0.3822504
  pull(Error) %>% 
  mean()

# Check mean values of relative bias in TFR
error_TFR_exp3B %>% 
 filter(Dataset == "25% Omission") %>% 19.75349
 # filter(Dataset == "100% Omission") %>% # 16.73654
  pull(Relative_Error) %>% 
  mean()

# Life Expectancy at birth ----
# Estimate life expectancy at birth from asmr 1x1 for the genealogical subsets ----
# with different proportions of omitted childless women

# Estimate age-specific mortality rates from the subset without 25% childless women
asmr_less_women_25_1 <- map_dfr(less_women_25, ~ estimate_mortality_rates_mod(opop = .x,
                                                                              final_sim_year = 2072, 
                                                                              year_min = 1750, 
                                                                              year_max = 2023, 
                                                                              year_group = 1,
                                                                              age_max_mort = 110, 
                                                                              age_group = 1),
                                .id = "Sim_id") 
save(asmr_less_women_25_1, file = "Measures/asmr_less_women_25_1.RData")

# Compute life tables from the subset without 25% childless women
lt_less_women_25 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_women_25_1)
save(lt_less_women_25, file = "Measures/lt_less_women_25.RData")


# Estimate age-specific mortality rates from the subset without 25% childless women
asmr_less_women_50_1 <- map_dfr(less_women_50, ~ estimate_mortality_rates_mod(opop = .x,
                                                                              final_sim_year = 2072,
                                                                              year_min = 1750, 
                                                                              year_max = 2023, 
                                                                              year_group = 1,
                                                                              age_max_mort = 110,
                                                                              age_group = 1),
                                .id = "Sim_id") 
save(asmr_less_women_50_1, file = "Measures/asmr_less_women_50_1.RData")

# Compute life tables from the subset without 50% childless women
lt_less_women_50 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_women_50_1)
save(lt_less_women_50, file = "Measures/lt_less_women_50.RData")


# Estimate age-specific mortality rates from the subset without 75% childless women
asmr_less_women_75_1 <- map_dfr(less_women_75, ~ estimate_mortality_rates_mod(opop = .x,
                                                                              final_sim_year = 2072, 
                                                                              year_min = 1750,
                                                                              year_max = 2023, 
                                                                              year_group = 1,
                                                                              age_max_mort = 110, 
                                                                              age_group = 1),
                                .id = "Sim_id") 
save(asmr_less_women_75_1, file = "Measures/asmr_less_women_75_1.RData")

# Compute life tables from the subset without 75% childless women
lt_less_women_75 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_women_75_1)
save(lt_less_women_75, file = "Measures/lt_less_women_75.RData")


# Estimate age-specific mortality rates from the subset without 100% childless women
asmr_less_women_100_1 <- map_dfr(less_women_100, ~ estimate_mortality_rates_mod(opop = .x,
                                                                                final_sim_year = 2072, 
                                                                                year_min = 1750,
                                                                                year_max = 2023, 
                                                                                year_group = 1,
                                                                                age_max_mort = 110, 
                                                                                age_group = 1),
                                 .id = "Sim_id") 
save(asmr_less_women_100_1, file = "Measures/asmr_less_women_100_1.RData")

# Compute life tables from the subset without 100% childless women
lt_less_women_100 <- lt_socsim_sims(asmr_socsim_sims = asmr_less_women_100_1)
save(lt_less_women_100, file = "Measures/lt_less_women_100.RData")

# Load and wrangle life tables for plotting ----

# Load life tables from each whole SOCSIM simulation
load("Measures/lt_10.RData")
# Load life tables from subset of All direct ancestors and their offspring
#load("Measures/lt_anc_off.RData")
load("Measures/lt_anc_off.RData")
# Load life tables from subset without 25% childless women
load("Measures/lt_less_women_25.RData")
# Load life tables from subset without 25% childless women
load("Measures/lt_less_women_50.RData")
# Load life tables from subset without 75% childless women
load("Measures/lt_less_women_75.RData")
# Load life tables from subset without 100% childless women
load("Measures/lt_less_women_100.RData")

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

# All direct ancestors and their offspring 
lt_anc_off2 <- lt_anc_off %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct Ancestors and their Offspring",
         Rate = "e0") %>%    
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All direct ancestors and their offspring without 25% childless women
lt_less_women_25b <- lt_less_women_25 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "25% Omission",
         Rate = "e0") %>%
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All direct ancestors and their offspring without 25% childless women
lt_less_women_50b <- lt_less_women_50 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "50% Omission",
         Rate = "e0") %>%    
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All direct ancestors and their offspring without 75% childless women
lt_less_women_75b <- lt_less_women_75 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "75% Omission",
         Rate = "e0") %>%    
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All direct ancestors and their offspring without 100% childless women
lt_less_women_100b <- lt_less_women_100 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "100% Omission",
         Rate = "e0") %>%    
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)


# Plot the estimates of life expectancy at birth
bind_rows(lt_whole2, lt_anc_off2,
          lt_less_women_25b, lt_less_women_50b, lt_less_women_75b, lt_less_women_100b) %>% 
  filter(Age == 0) %>%
  mutate(Dataset = factor(Dataset, levels = c("100% Omission", "75% Omission", "50% Omission", "25% Omission", 
                                              "Direct Ancestors and their Offspring", "Whole Simulation"))) %>% 
  group_by(Year, sex, Dataset) %>% 
  summarise(ex = mean(ex, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year, y = ex, group = Dataset, colour = Dataset)) +
  facet_grid(. ~ sex) +
  geom_line(linewidth = 1.3)+ 
  scale_color_viridis(option = "D", discrete = T, direction = -1)+
  theme_graphs() 
ggsave(file="Graphs/Socsim_Exp3B_e0.jpeg", width=17, height=9, dpi=300)

#----------------------------------------------------------------------------------------------------
# Figure for EPC presentation

# Plot the estimates of life expectancy at birth
bind_rows(lt_whole2,  lt_anc_off2, lt_less_women_25b, lt_less_women_100b)  %>%
  filter(Age == 0) %>%
  mutate(Sex = ifelse(sex == "female", "Female", "Male")) %>%
  group_by(Year, Sex, Dataset) %>% 
  summarise(ex = mean(ex, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year, y = ex, group = Dataset))+
  geom_line(aes(colour = Dataset), linewidth = 1.3)+
  scale_color_manual(values = c( "#FF834C","#E7495B","#75007A", "#007A75"))+
  facet_wrap(~Sex) +
  theme_graphs() +
  labs(y = "Life expectancy at birth")

ggsave(file="Graphs/Socsim_Exp3B_e0_grp.jpeg", width=17, height=9, dpi=300)
#----------------------------------------------------------------------------------------------------
# Summary measure of error in e0 ----

# Difference in means
DiM_e0_Exp3B <- bind_rows(lt_whole2, 
                   lt_less_women_25b, lt_less_women_50b, lt_less_women_75b, lt_less_women_100b) %>%
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
MoD_e0_Exp3B <- bind_rows(lt_whole2, 
                   lt_less_women_25b, lt_less_women_50b, lt_less_women_75b, lt_less_women_100b) %>%
  filter(Year > 1750 & Age == 0) %>% 
  pivot_wider(id_cols = c(Year, sex, Sim_id), names_from = "Dataset", values_from = "ex") %>% 
  pivot_longer(cols = 5:8, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation`, 
         Relative_Error = (Error/`Whole Simulation`)*100) %>% 
  group_by(Year, sex, Dataset) %>% 
  summarise(Error = mean(Error, na.rm = T), 
          Relative_Error = mean(Relative_Error, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Type = "MoD")

# Bind both error measures and save the data frame
error_e0_exp3B <- bind_rows(DiM_e0_Exp3B, MoD_e0_Exp3B) 
save(error_e0_exp3B, file = "Measures/error_e0_exp3B.RData")

# Absolute error
error_e0_exp3B %>% 
  ggplot(aes(x = Year, y = Error, group = Dataset, colour = Type)) +
  facet_grid(Dataset ~ sex)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3B_e0_Error.jpeg", width=17, height=9, dpi=300)

# Relative error
error_e0_exp3B %>% 
  ggplot(aes(x = Year, y = Relative_Error, group = Dataset, colour = Type)) +
  facet_grid(Dataset ~ sex)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3B_e0_Rel_Error.jpeg", width=17, height=9, dpi=300)


# Check mean values of bias in e0 over the whole period
error_e0_exp3B %>% 
  filter(sex == "female") %>%
  #filter(Dataset == "25% Omission") %>% # 1.050917
  filter(Dataset == "100% Omission") %>% # 1.493743
  pull(Error) %>%
  mean()


# Check mean values of relative bias in e0 over the whole period for 100% omission
error_e0_exp3B %>% 
  filter(sex == "female") %>%
 # filter(Dataset == "25% Omission") %>% # 2.113992
  filter(Dataset == "100% Omission") %>% # 3.059821
  pull(Relative_Error) %>%
  mean()

#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

yrs_plot2 <- c(1750, 1800, 1850, 1900, 1950, 2000)

# Define the same y breaks for all plots
y_breaks_TFR <- c(0:5)
y_breaks_e0 <- c(20, 40, 60, 80)

## TFR and e0 (for females) from whole SOCSIM simulation and subsets of direct ancestors and their offspring

Summary_Exp3B <- 
bind_rows(TFR_whole %>% rename(Estimate = TFR),
          TFR_anc_off %>% rename(Estimate = TFR),
          TFR_less_women_25 %>% rename(Estimate = TFR),
          # TFR_less_women_50 %>% rename(Estimate = TFR),
          # TFR_less_women_75 %>% rename(Estimate = TFR),
          TFR_less_women_100 %>% rename(Estimate = TFR)) %>%
  mutate(sex = "female") %>%  
  bind_rows(lt_whole2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_off2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_less_women_25b %>% rename(Estimate = ex) %>% filter(Age == 0),
            # lt_less_women_50b %>% rename(Estimate = ex) %>% filter(Age == 0),
            # lt_less_women_75b %>% rename(Estimate = ex) %>% filter(Age == 0), 
            lt_less_women_100b %>% rename(Estimate = ex) %>% filter(Age == 0)) %>%
  filter(sex == "female") %>%
  group_by(Year, Dataset, Rate, sex) %>% 
  summarise(Estimate = mean(Estimate, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life Expectancy at Birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life Expectancy at Birth")),
         Dataset = factor(Dataset, levels = c("100% Omission", "25% Omission", 
                                              "Direct Ancestors and their Offspring", "Whole Simulation"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2),
             aes(shape = Dataset), size = 11) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c( "#FF834C","#E7495B","#75007A", "#007A75"))+
  scale_shape_manual(values = c(8, 22, 18, 46)) + 
  scale_x_continuous(breaks = yrs_plot2)+
  facetted_pos_scales(y = list("Total Fertility Rate" = scale_y_continuous(breaks = y_breaks_TFR, 
                                                                           limits = c(0, NA)),
                               "Life Expectancy at Birth" =  scale_y_continuous(breaks = y_breaks_e0)))+
  theme_graphs() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
Summary_Exp3B
# Save the plot
ggsave(file="Graphs/Socsim_Exp3B_TFR_e0.jpeg", width=17, height=9, dpi=300)

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

By_Age_Exp3B + 
  geom_text(data = plot_labs1, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") + 
  theme(plot.margin = margin(0,0,1,0, "cm")) +
  Summary_Exp3B + 
  geom_text(data = plot_labs2, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") +
  plot_layout(ncol = 1)
ggsave(file="Final_Graphs/Final_Socsim_Exp3B_Combined.jpeg", width=18, height=21, dpi=300)