#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Sweden - Trace and compare subsets with randomly removed proportions of childless women
# U:/SOCSIM/SOCSIM_Genealogies/6_Compare_Omitted_women.R

# Remove a proportion of sub-populations (childless women)
# from SOCSIM microsimulations for Sweden (1751-2022) 
# Trace genealogies and compare demographic measures from the whole simulation and the subsets 

# Created by Liliana Calderon on 11-07-2023
# Last modified by Liliana Calderon on 02-08-2023

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
source("Functions_Graphs.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with asmr calculated with rsocsim::estimate_mortality_rates()
source("Functions_Life_Table.R")

#------------------------------------------------------------------------------------------------------
## Load necessary data and randomly removed proportions of childless women ----

# Load the data frame with the ancestors and collateral kin of 10 simulations samples
load("Subsets/anc_kin_10.RData")

# NB: This data still contains duplicates within each simulation, as someone can be a relative of different egos
# As we do not distinguish now between type of kin, we will keep only unique pids before merging. 
anc_kin_10 <- anc_kin_10 %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup()  
anc_kin_list <- anc_kin_10 %>% split(.$Sim_id)


# Function to get a sample of childless women, born after 1735
# who survived at least until reproductive ages

sample_women <- function(opop = opop, age_threshold, percentage) {
  
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
#  Estimate ASFR and ASMR for the genealogical subsets with direct ancestors and collateral kin, 
# after removing a proportion of childless women

# All direct ancestors and kin from the subset without 25% childless women
# Create a list of data frames with opop of the filter data
less_women_25 <- less_women_25 %>% split(.$Sim_id)

# Estimate age-specific fertility rates from the subset without 25% childless women
asfr_less_women_25 <- map_dfr(less_women_25, ~ estimate_fertility_rates(opop = .x,
                                                                        final_sim_year = 2022, 
                                                                        year_min = 1750, 
                                                                        year_max = 2020, 
                                                                        year_group = 5,
                                                                        age_min_fert = 10, 
                                                                        age_max_fert = 55, 
                                                                        age_group = 5), 
                              .id = "Sim_id") 
save(asfr_less_women_25, file = "Measures/asfr_less_women_25.RData")

# Estimate age-specific mortality rates from the subset without 25% childless women
asmr_less_women_25 <- map_dfr(less_women_25, ~ estimate_mortality_rates(opop = .x,
                                                                        final_sim_year = 2022, 
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
asfr_less_women_50 <- map_dfr(less_women_50, ~ estimate_fertility_rates(opop = .x,
                                                                        final_sim_year = 2022,
                                                                        year_min = 1750, 
                                                                        year_max = 2020, 
                                                                        year_group = 5,
                                                                        age_min_fert = 10, 
                                                                        age_max_fert = 55, 
                                                                        age_group = 5), 
                              .id = "Sim_id") 
save(asfr_less_women_50, file = "Measures/asfr_less_women_50.RData")

# Estimate age-specific mortality rates from the subset without 25% childless women
asmr_less_women_50 <- map_dfr(less_women_50, ~ estimate_mortality_rates(opop = .x,
                                                                        final_sim_year = 2022,
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
asfr_less_women_75 <- map_dfr(less_women_75, ~ estimate_fertility_rates(opop = .x,
                                                                        final_sim_year = 2022, 
                                                                        year_min = 1750, 
                                                                        year_max = 2020,
                                                                        year_group = 5,
                                                                        age_min_fert = 10, 
                                                                        age_max_fert = 55, 
                                                                        age_group = 5), 
                              .id = "Sim_id") 
save(asfr_less_women_75, file = "Measures/asfr_less_women_75.RData")

# Estimate age-specific mortality rates from the subset without 75% childless women
asmr_less_women_75 <- map_dfr(less_women_75, ~ estimate_mortality_rates(opop = .x,
                                                                        final_sim_year = 2022, 
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
asfr_less_women_100 <- map_dfr(less_women_100, ~ estimate_fertility_rates(opop = .x,
                                                                          final_sim_year = 2022, 
                                                                          year_min = 1750, 
                                                                          year_max = 2020,
                                                                          year_group = 5,
                                                                          age_min_fert = 10, 
                                                                          age_max_fert = 55, 
                                                                          age_group = 5), 
                               .id = "Sim_id") 
save(asfr_less_women_100, file = "Measures/asfr_less_women_100.RData")

# Estimate age-specific mortality rates from the subset without 100% childless women
asmr_less_women_100 <- map_dfr(less_women_75, ~ estimate_mortality_rates(opop = .x,
                                                                         final_sim_year = 2022, 
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
# Load asfr from the subset of all direct ancestors and collateral kin
load("Measures/asfr_anc_col.RData")
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

# All Direct ancestors and collateral kin 
asfr_anc_col2 <- asfr_anc_col %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "ASFR") 

# All Direct ancestors and collateral kin without 25% childless women
asfr_less_women_25b <- asfr_less_women_25 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "25% Omission",
         Rate = "ASFR")

# All Direct ancestors and collateral kin without 50% childless women
asfr_less_women_50b <- asfr_less_women_50 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "50% Omission",
         Rate = "ASFR")

# All Direct ancestors and collateral kin without 75% childless women
asfr_less_women_75b <- asfr_less_women_75 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "75% Omission",
         Rate = "ASFR")

# All Direct ancestors and collateral kin without 100% childless women
asfr_less_women_100b <- asfr_less_women_100 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "100% Omission",
         Rate = "ASFR")

## Plot ASFR from whole simulation with genealogical subsets with different proportions of missing childless women

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2, asfr_less_women_25b, asfr_less_women_50b, asfr_less_women_75b, asfr_less_women_100b) %>% 
  filter(year %in% yrs_plot & !is.nan(ASFR)) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset), colour = year))+
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(shape = Dataset), size = 6)+ 
  scale_color_manual(values = c("#79B727","#B72779", "#2779B7"))+ 
  scale_shape_manual(values = c(8, 21, 22, 23, 46)) +
  theme_graphs()
# labs(title = "Age-Specific Fertility Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and genealogical subsets with missing women) 
ggsave(file="Graphs/Socsim_Exp3B_ASFR.jpeg", width=17, height=9, dpi=200)


# ASMR ----

# Load ASMR 5x5 from the 10 simulations
load("Measures/asmr_10.RData")
# Load asmr from the subset of all direct ancestors and collateral kin, calculated on 3_Compare_Ancestors
load("Measures/asmr_anc_col.RData")
# Load asmr for the genealogical subset without 25% childless women
load("Measures/asmr_less_women_25.RData")
# Load asmr for the genealogical subset without 25% childless women
load("Measures/asmr_less_women_50.RData")
# Load asmr for the genealogical subset without 75% childless women
load("Measures/asmr_less_women_75.RData")
# Load asmr for the genealogical subset without 100% childless women
load("Measures/asmr_less_women_100.RData")

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_anc_col$age)

# Whole SOCSIM simulations
asmr_whole2 <- asmr_10 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "Whole Simulation", 
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

asmr_anc_col2 <- asmr_anc_col %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "ASMR", 
         Omitted = "NA") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All Direct ancestors and collateral kin without 25% childless women
asmr_less_women_25b <- asmr_less_women_25 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "25% Omission",
         Rate = "ASMR") %>%
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All Direct ancestors and collateral kin without 50% childless women
asmr_less_women_50b <- asmr_less_women_50 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "50% Omission",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All Direct ancestors and collateral kin without 75% childless women
asmr_less_women_75b <- asmr_less_women_75 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "75% Omission",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All Direct ancestors and collateral kin without 75% childless women
asmr_less_women_100b <- asmr_less_women_100 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "100% Omission",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)


## Plotting ASMR from whole simulation with genealogical subsets with different proportions of missing childless women

# Same years to plot than above (in intervals). Change if necessary
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
  scale_color_manual(values = c("#79B727","#B72779", "#2779B7"))+ 
  scale_shape_manual(values = c(8, 21, 22, 23, 46)) +
  theme_graphs()+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs()
#labs(title = "Age-Specific Mortality Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and genealogical subsets with omitted women") 
ggsave(file="Graphs/Socsim_Exp3B_ASMR.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Final plot combining ASFR and ASMR ----

# Years to plot
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_whole2$age)

## Plotting ASFR and ASMR (for females) from whole SOCSIM simulation and genealogical subsets
By_Age_Exp3B <- 
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_anc_col2 %>% rename(Estimate = ASFR), 
          asfr_less_women_25b %>% rename(Estimate = ASFR),
          asfr_less_women_50b %>% rename(Estimate = ASFR),
          asfr_less_women_75b %>% rename(Estimate = ASFR), 
          asfr_less_women_100b %>% rename(Estimate = ASFR)) %>%
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx), 
            asmr_anc_col2 %>% rename(Estimate = mx), 
            asmr_less_women_25b %>% rename(Estimate = mx),
            asmr_less_women_50b %>% rename(Estimate = mx),
            asmr_less_women_75b %>% rename(Estimate = mx),
            asmr_less_women_100b %>% rename(Estimate = mx)) %>%
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
  scale_shape_manual(values = c(8, 23,  22, 21, 18, 46)) + 
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() + 
  labs(x = "Age") +
  guides(shape = guide_legend(order = 1), col = guide_legend(order = 2)) +
  theme(legend.justification = "left", 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
By_Age_Exp3B
ggsave(file="Graphs/Final_Socsim_Exp3B_ASFR_ASMR.jpeg", width=17, height=9, dpi=200)
#----------------------------------------------------------------------------------------------------
## Summary measures: TFR and e0 ----
# Here, we use the rates by 1 year age group and 1 calendar year

# Total Fertility Rate ----
# Estimate Total Fertility Rate from asfr 1x1 ----

# Estimate age-specific fertility rates 1x1 from the subset without 25% childless women
asfr_less_women_25_1 <- map_dfr(less_women_25, ~ estimate_fertility_rates(opop = .x,
                                                                          final_sim_year = 2022,
                                                                          year_min = 1750,
                                                                          year_max = 2023,
                                                                          year_group = 1,
                                                                          age_min_fert = 10,
                                                                          age_max_fert = 55, 
                                                                          age_group = 1), 
                                .id = "Sim_id") 
save(asfr_less_women_25_1, file = "Measures/asfr_less_women_25_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 25% childless women
asfr_less_women_50_1 <- map_dfr(less_women_50, ~ estimate_fertility_rates(opop = .x,
                                                                          final_sim_year = 2022,
                                                                          year_min = 1750, 
                                                                          year_max = 2023, 
                                                                          year_group = 1,
                                                                          age_min_fert = 10, 
                                                                          age_max_fert = 55, 
                                                                          age_group = 1), 
                                .id = "Sim_id") 
save(asfr_less_women_50_1, file = "Measures/asfr_less_women_50_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 75% childless women
asfr_less_women_75_1 <- map_dfr(less_women_75, ~ estimate_fertility_rates(opop = .x,
                                                                          final_sim_year = 2022, 
                                                                          year_min = 1750, 
                                                                          year_max = 2023, 
                                                                          year_group = 1,
                                                                          age_min_fert = 10, 
                                                                          age_max_fert = 55, 
                                                                          age_group = 1), 
                                .id = "Sim_id") 
save(asfr_less_women_75_1, file = "Measures/asfr_less_women_75_1.RData")

# Estimate age-specific fertility rates 1x1 from the subset without 100% childless women
asfr_less_women_100_1 <- map_dfr(less_women_100, ~ estimate_fertility_rates(opop = .x,
                                                                            final_sim_year = 2022, 
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
# Load asfr 1x1 from the subset with all direct ancestors and collateral kin
load("Measures/asfr_anc_col_1.RData")
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

# All Direct ancestors and collateral kin without 25% childless women
TFR_less_women_25 <- asfr_less_women_25_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "25% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 50% childless women
TFR_less_women_50 <- asfr_less_women_50_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "50% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 75% childless women
TFR_less_women_75 <- asfr_less_women_75_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "75% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 100% childless women
TFR_less_women_100 <- asfr_less_women_100_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "100% Omission",
         Rate = "TFR",           
         sex = "female")

## Plot TFR from whole SOCSIM simulation and subsets with different proportions of omitted children

bind_rows(TFR_whole, TFR_anc_col,
          TFR_less_women_25, TFR_less_women_50, TFR_less_women_75, TFR_less_women_100) %>% 
  group_by(Year, Dataset) %>% 
  summarise(TFR = mean(TFR, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year, y = TFR, group = Dataset, colour = Dataset)) +
  geom_line(linewidth = 1.3)+ 
  scale_color_viridis(option = "D", discrete = T, direction = -1)+
  theme_graphs() 
# labs(title = "Total Fertility Rate in Sweden (1751-2022), 
#retrieved from a SOCSIM simulation and subsets with different proportions of omitted women") 
ggsave(file="Graphs/Socsim_Exp3B_TFR.jpeg", width=17, height=9, dpi=200)


# Summary measure of error in TFR ----

# Differences of means
DM_TFR_Exp3B <- bind_rows(TFR_whole, TFR_anc_col,
           TFR_less_women_25, TFR_less_women_50, TFR_less_women_75, TFR_less_women_100) %>%
  filter(Year > 1750) %>% 
  group_by(Year, Dataset) %>% 
  summarise(TFR = mean(TFR, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = Year, names_from = "Dataset", values_from = "TFR") %>% 
  pivot_longer(cols = 2:6, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = `Whole Simulation` - Genealogy, 
         Type = "DM")

# Mean of differences
MD_TFR_Exp3B <- bind_rows(TFR_whole, TFR_anc_col,
                    TFR_less_women_25, TFR_less_women_50, TFR_less_women_75, TFR_less_women_100) %>%
  filter(Year > 1750) %>% 
  pivot_wider(id_cols = c(Year, Sim_id), names_from = "Dataset", values_from = "TFR") %>% 
  pivot_longer(cols = 4:8, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = `Whole Simulation` - Genealogy) %>% 
  group_by(Year, Dataset) %>% 
  summarise(Error = mean(Error, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Type = "MD")

bind_rows(DM_TFR_Exp3B, MD_TFR_Exp3B) %>%
  ggplot(aes(x = Year, y = Error, colour = Dataset, group = Dataset)) +
  geom_line(linewidth = 1.3)+
  geom_point(aes(shape = Type), size = 3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3B_TFR_Error.jpeg", width=17, height=9, dpi=200)

# Life Expectancy at birth ----
# Estimate life expectancy at birth from asmr 1x1 for the genealogical subsets ----
# with different proportions of omitted childless women

# Estimate age-specific mortality rates from the subset without 25% childless women
asmr_less_women_25_1 <- map_dfr(less_women_25, ~ estimate_mortality_rates(opop = .x,
                                                                          final_sim_year = 2022, 
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
asmr_less_women_50_1 <- map_dfr(less_women_50, ~ estimate_mortality_rates(opop = .x,
                                                                          final_sim_year = 2022,
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
asmr_less_women_75_1 <- map_dfr(less_women_75, ~ estimate_mortality_rates(opop = .x,
                                                                          final_sim_year = 2022, 
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
asmr_less_women_100_1 <- map_dfr(less_women_100, ~ estimate_mortality_rates(opop = .x,
                                                                            final_sim_year = 2022, 
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
# Load life tables from subset of all direct ancestors and collateral kin
load("Measures/lt_anc_col.RData")
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

# All Direct ancestors and collateral kin 
lt_anc_col2 <- lt_anc_col %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "e0") %>%    
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 25% childless women
lt_less_women_25b <- lt_less_women_25 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "25% Omission",
         Rate = "e0") %>%
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 25% childless women
lt_less_women_50b <- lt_less_women_50 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "50% Omission",
         Rate = "e0") %>%    
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 75% childless women
lt_less_women_75b <- lt_less_women_75 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "75% Omission",
         Rate = "e0") %>%    
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 100% childless women
lt_less_women_100b <- lt_less_women_100 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "100% Omission",
         Rate = "e0") %>%    
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)


# Plot the estimates of life expectancy at birth
bind_rows(lt_whole2, lt_anc_col2,
          lt_less_women_25b, lt_less_women_50b, lt_less_women_75b, lt_less_women_100b) %>% 
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
# labs(title = ""Life expectancy at birth in Sweden (e0), 1751-2020,
#retrieved from a SOCSIM simulation and subsets with different proportions of omitted women") 
ggsave(file="Graphs/Socsim_Exp3B_e0.jpeg", width=17, height=9, dpi=200)


# Summary measure of error in e0 ----

# Differences of means
DM_e0_Exp3B <- bind_rows(lt_whole2, 
                   lt_less_women_25b, lt_less_women_50b, lt_less_women_75b, lt_less_women_100b) %>%
  filter(Year > 1750 & Age == 0) %>% 
  group_by(Year, sex, Dataset) %>% 
  summarise(ex = mean(ex, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(Year:sex), names_from = "Dataset", values_from = "ex") %>% 
  pivot_longer(cols = 3:6, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = `Whole Simulation` - Genealogy, 
         Type = "DM") %>% 
  select(Year, sex, Dataset, Error, Type) 

# Mean of differences
MD_e0_Exp3B <- bind_rows(lt_whole2, 
                   lt_less_women_25b, lt_less_women_50b, lt_less_women_75b, lt_less_women_100b) %>%
  filter(Year > 1750 & Age == 0) %>% 
  pivot_wider(id_cols = c(Year, sex, Sim_id), names_from = "Dataset", values_from = "ex") %>% 
  pivot_longer(cols = 5:8, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = `Whole Simulation` - Genealogy) %>% 
  group_by(Year, sex, Dataset) %>% 
  summarise(Error = mean(Error, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Type = "MD")

bind_rows(DM_e0_Exp3B, MD_e0_Exp3B) %>%
  ggplot(aes(x = Year, y = Error, colour = Dataset, group = Dataset)) +
  facet_wrap(. ~ sex)+
  geom_line(linewidth = 1.3)+
  geom_point(aes(shape = Type), size = 3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp3B_e0_Error.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

yrs_plot2 <- c(1750, 1800, 1850, 1900, 1950, 2000)

## TFR and e0 (for females) from whole SOCSIM simulation and subsets of "direct" and all collateral kin

Summary_Exp3B <- 
bind_rows(TFR_whole %>% rename(Estimate = TFR),
          TFR_anc_col %>% rename(Estimate = TFR),
          TFR_less_women_25 %>% rename(Estimate = TFR),
          TFR_less_women_50 %>% rename(Estimate = TFR),
          TFR_less_women_75 %>% rename(Estimate = TFR),
          TFR_less_women_100 %>% rename(Estimate = TFR)) %>%
  mutate(sex = "female") %>%  
  bind_rows(lt_whole2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_col2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_less_women_25b %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_less_women_50b %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_less_women_75b %>% rename(Estimate = ex) %>% filter(Age == 0), 
            lt_less_women_100b %>% rename(Estimate = ex) %>% filter(Age == 0)) %>%
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
  scale_color_manual(values = c("#FEBD2A", "#F48849", "#DB5C68", "#B83289",  "#75007A", "#007A75"))+
  scale_shape_manual(values = c(8, 23,  22, 21, 18, 46)) +
  theme_graphs() +
  theme(legend.justification = "left",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
Summary_Exp3B
# labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and subsets with different proportions of omitted children")
# Save the plot
ggsave(file="Graphs/Final_Socsim_Exp3B_TFR_e0.jpeg", width=17, height=9, dpi=200)

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
ggsave(file="Graphs/Final_Socsim_Exp3B_Combined.jpeg", width=18, height=21, dpi=200)