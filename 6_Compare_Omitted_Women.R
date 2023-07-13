#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Sweden - Trace and compare subsets with randomly removed proportions of childless women
# U:/SOCSIM/SOCSIM_Genealogies/6_Compare_Omitted_women.R

# Remove a proportion of sub-populations (childless women)
# from SOCSIM microsimulations for Sweden (1751-2022) 
# Trace genealogies and compare demographic measures from the whole simulation and the subsets 

# Created by Liliana Calderon on 11-07-2023
# Last modified by Liliana Calderon on 12-07-2023

## NB: To run this code, it is necessary to have already run the scripts 
# 1_Run_Simulations.R, 3_Compare_Ancestors.R and 4_Compare_Kin.R
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

#------------------------------------------------------------------------------------------------------
## Load data with simulation results, sample of egos alive in 2023, subset of ancestors and lateral kin 

# Load saved list with opop from 10 simulations, generated in 1_Run_Simulations.R
load("sims_opop.RData")

# Load the data frame with the ancestors and relatives of 10 simulations samples
load("Subsets/anc_kin_10.RData")

# NB: This data still contains duplicates within each simulation, as someone can be a relative of different egos
# As we do not distinguish now between type of kin, we will keep only unique pids before merging. 
# We remove the grand-children as they are not included in the genealogies
anc_kin_10 <- anc_kin_10 %>% 
  filter(!kin_type %in% c("gchildren")) %>%
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup()


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

# Get samples of childless women for each whole simulation, using different proportions of omission
# and remove them from the anc_kin opop
# NB: Not all these childless women are included in a genealogy. 
# Hence, population in less_women is not equal to anc_kin_10 - miss_women

miss_women_05 <- map_dfr(sims_opop, ~ sample_women(opop = .x, percentage = 5),
                             .id = "Sim_id") 
less_women_05 <- anti_join(anc_kin_10, miss_women_05)
save(less_women_05, file = "Subsets/less_women_05.RData")

miss_women_10 <- map_dfr(sims_opop, ~ sample_women(opop = .x, percentage = 10),
                          .id = "Sim_id")
less_women_10 <- anti_join(anc_kin_10, miss_women_10)
save(less_women_10, file = "Subsets/less_women_10.RData")

miss_women_20 <- map_dfr(sims_opop, ~ sample_women(opop = .x, percentage = 20),
                          .id = "Sim_id") 
less_women_20 <- anti_join(anc_kin_10, miss_women_20) 
save(less_women_20, file = "Subsets/less_women_20.RData")

#----------------------------------------------------------------------------------------------------
## Age-Specific Fertility and Mortality rates, 5x5  -----
#  Estimate ASFR and ASMR for the genealogical subsets with direct ancestors and collateral kin, 
# after removing a proportion of childless women

# All direct ancestors and kin from the subset without 5% childless women
# Create a list of data frames with opop of the filter data
less_women_05 <- less_women_05 %>% split(.$Sim_id)

# Estimate age-specific fertility rates for the subset without 5% childless women
asfr_less_women_05 <- map_dfr(less_women_05, ~ estimate_fertility_rates(opop = .x,
                                                                        final_sim_year = 2022, 
                                                                        year_min = 1750, 
                                                                        year_max = 2020, 
                                                                        year_group = 5,
                                                                        age_min_fert = 10, 
                                                                        age_max_fert = 55, 
                                                                        age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_women_05, file = "Measures/asfr_less_women_05.RData")

# Estimate age-specific mortality rates for the subset without 5% childless women
asmr_less_women_05 <- map_dfr(less_women_05, ~ estimate_mortality_rates(opop = .x,
                                                                        final_sim_year = 2022, 
                                                                        year_min = 1750, 
                                                                        year_max = 2020, 
                                                                        year_group = 5,
                                                                        age_max_mort = 110,
                                                                        age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_women_05, file = "Measures/asmr_less_women_05.RData")


# All direct ancestors and kin from the subset without 10% childless women
# Create a list of data frames with opop of the filter data
less_women_10 <- less_women_10 %>% split(.$Sim_id)

# Estimate age-specific fertility rates for the subset without 10% childless women
asfr_less_women_10 <- map_dfr(less_women_10, ~ estimate_fertility_rates(opop = .x,
                                                                        final_sim_year = 2022,
                                                                        year_min = 1750, 
                                                                        year_max = 2020, 
                                                                        year_group = 5,
                                                                        age_min_fert = 10, 
                                                                        age_max_fert = 55, 
                                                                        age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_women_10, file = "Measures/asfr_less_women_10.RData")

# Estimate age-specific mortality rates for the subset without 10% childless women
asmr_less_women_10 <- map_dfr(less_women_10, ~ estimate_mortality_rates(opop = .x,
                                                                        final_sim_year = 2022,
                                                                        year_min = 1750, 
                                                                        year_max = 2020, 
                                                                        year_group = 5,
                                                                        age_max_mort = 110, 
                                                                        age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_women_10, file = "Measures/asmr_less_women_10.RData")


# All direct ancestors and kin from the subset without 20% childless women
# Create a list of data frames with opop of the filter data
less_women_20 <- less_women_20 %>% split(.$Sim_id)

# Estimate age-specific fertility rates for the subset without 20% childless women
asfr_less_women_20 <- map_dfr(less_women_20, ~ estimate_fertility_rates(opop = .x,
                                                                        final_sim_year = 2022, 
                                                                        year_min = 1750, 
                                                                        year_max = 2020,
                                                                        year_group = 5,
                                                                        age_min_fert = 10, 
                                                                        age_max_fert = 55, 
                                                                        age_group = 5), 
                                   .id = "Sim_id") 
save(asfr_less_women_20, file = "Measures/asfr_less_women_20.RData")

# Estimate age-specific mortality rates for the subset without 20% childless women
asmr_less_women_20 <- map_dfr(less_women_20, ~ estimate_mortality_rates(opop = .x,
                                                                        final_sim_year = 2022, 
                                                                        year_min = 1750, 
                                                                        year_max = 2020, 
                                                                        year_group = 5,
                                                                        age_max_mort = 110, 
                                                                        age_group = 5),
                                   .id = "Sim_id") 
save(asmr_less_women_20, file = "Measures/asmr_less_women_20.RData")


#----------------------------------------------------------------------------------------------------
## Plot estimates from genealogical subsets with some children removed ----

# Load asfr for the subset of all direct ancestors and collateral kin
load("Measures/asfr_anc_zaukgausc.RData")
# Load asmr for the subset of all direct ancestors and collateral kin
load("Measures/asmr_anc_zaukgausc.RData")

# Load asfr for the genealogical subset without 5% childless women
load("Measures/asfr_less_women_05.RData")
# Load asmr for the genealogical subset without 5% childless women
load("Measures/asmr_less_women_05.RData")
# Load asfr for the genealogical subset without 10% childless women
load("Measures/asfr_less_women_10.RData")
# Load asmr for the genealogical subset without 10% childless women
load("Measures/asmr_less_women_10.RData")
# Load asfr for the genealogical subset without 20% childless women
load("Measures/asfr_less_women_20.RData")
# Load asmr for the genealogical subset without 20% childless women
load("Measures/asmr_less_women_20.RData")


# Choose years to plot (in intervals).
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_anc_zaukgausc$age)

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

## Plot ASFR and ASMR (for women) for the subset with all direct ancestors and collateral kin
plot_asfr_asmr(asfr_anc_zaukgausc, asmr_anc_zaukgausc, yrs_plot, age_levels)
# Plot ASFR and ASMR (for women) for the subset without 5% childless women
plot_asfr_asmr(asfr_less_women_05, asmr_less_women_05, yrs_plot, age_levels)
# Plot ASFR and ASMR (for women) for the subset without 10% childless women
plot_asfr_asmr(asfr_less_women_10, asmr_less_women_10, yrs_plot, age_levels)
# Plot ASFR and ASMR (for women) for the subset without 20% childless women
plot_asfr_asmr(asfr_less_women_20, asmr_less_women_20, yrs_plot, age_levels)

#----------------------------------------------------------------------------------------------------
## Comparison of whole simulation with genealogical subsets with different proportions of missing childless women ----

# ASFR ----

# Load mean ASFR 5x5 from the 10 simulations, calculated on 2_Compare_Input_Output
load("Measures/asfr_whole.RData")

## Calculate the mean of the different simulations and add relevant columns

# Whole SOCSIM simulation
asfr_whole2 <- asfr_whole %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Whole Simulation", 
         Rate = "ASFR") 

# All Direct ancestors and collateral kin 
# Not included
asfr_anc_zaukgausc2 <- asfr_anc_zaukgausc %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "ASFR") 

# All Direct ancestors and collateral kin without 5% childless women
asfr_less_women_05b <- asfr_less_women_05 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "05% Omission",
         Rate = "ASFR")

# All Direct ancestors and collateral kin without 10% childless women
asfr_less_women_10b <- asfr_less_women_10 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "10% Omission",
         Rate = "ASFR")

# All Direct ancestors and collateral kin without 20% childless women
asfr_less_women_20b <- asfr_less_women_20 %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "20% Omission",
         Rate = "ASFR")

## Plot ASFR from whole simulation with genealogical subsets with different proportions of missing childless women

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2, asfr_less_women_05b, asfr_less_women_10b, asfr_less_women_20b) %>% 
  filter(year %in% yrs_plot & !is.nan(ASFR)) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset), colour = year))+
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(shape = Dataset), size = 6)+ 
  scale_color_manual(values = c("#B72779", "#2779B7"))+
  scale_shape_manual(values = c(21, 22, 23, 46)) +
  theme_graphs()
# labs(title = "Age-Specific Fertility Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and genealogical subsets with missing women) 
ggsave(file="Graphs/Socsim_Exp4_ASFR.jpeg", width=17, height=9, dpi=200)


# ASMR ----

# Load mean ASMR rates 5x5 from the 10 simulations, calculated on 2_Compare_Input_Output
load("Measures/asmr_whole.RData")

# Whole SOCSIM simulation
asmr_whole2 <- asmr_whole %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "Whole Simulation", 
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All Direct ancestors and collateral kin 
asmr_anc_zaukgausc2 <- asmr_anc_zaukgausc %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "ASMR", 
         Omitted = "NA") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All Direct ancestors and collateral kin without 5% childless women
asmr_less_women_05b <- asmr_less_women_05 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "05% Omission",
         Rate = "ASMR") %>%
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All Direct ancestors and collateral kin without 10% childless women
asmr_less_women_10b <- asmr_less_women_10 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "10% Omission",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# All Direct ancestors and collateral kin without 20% childless women
asmr_less_women_20b <- asmr_less_women_20 %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),          
         Dataset = "20% Omission",
         Rate = "ASMR") %>%    
  select(year, age, Sex, mx = socsim, Dataset, Rate)


## Plotting ASMR from whole simulation with genealogical subsets with different proportions of missing childless women

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1900,1905)", "[2000,2005)") 

bind_rows(asmr_whole2, asmr_less_women_05b, asmr_less_women_10b, asmr_less_women_20b) %>% 
  rename(Year = year) %>% 
  filter(Year %in% yrs_plot) %>%
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>% 
  ggplot(aes(x = age, y = mx, group = interaction(Year, Dataset), colour = Year))+
  facet_grid(. ~ Sex) +
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(shape = Dataset), size = 11)+ 
  scale_color_manual(values = c("#B72779", "#2779B7"))+
  scale_shape_manual(values = c(21, 22, 23, 46)) +
  theme_graphs()+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs()
#labs(title = "Age-Specific Mortality Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and genealogical subsets with omitted women") 
ggsave(file="Graphs/Socsim_Exp4_ASMR.jpeg", width=17, height=9, dpi=200)


#----------------------------------------------------------------------------------------------------
## Final plot combining ASFR and ASMR ----

# Years to plot limited to  two years
yrs_plot <- c("[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_whole2$age)

## Plotting ASFR and ASMR (for females) from whole SOCSIM simulation and genealogical subsets
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_anc_zaukgausc2 %>% rename(Estimate = ASFR), 
          asfr_less_women_05b %>% rename(Estimate = ASFR),
          asfr_less_women_10b %>% rename(Estimate = ASFR),
          asfr_less_women_20b %>% rename(Estimate = ASFR)) %>%
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx), 
            asmr_anc_zaukgausc2 %>% rename(Estimate = mx), 
            asmr_less_women_05b %>% rename(Estimate = mx),
            asmr_less_women_10b %>% rename(Estimate = mx),
            asmr_less_women_20b %>% rename(Estimate = mx)) %>%
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels),
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates")) %>%
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(linewidth = 1.2, show.legend = T)+ 
  geom_point(aes(shape = Dataset), size = 11)+ 
  scale_color_manual(values = c("#B72779", "#2779B7"))+
  scale_shape_manual(values = c(21, 22, 23, 18, 46)) + 
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() + 
  theme(legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15))+
  labs(x = "Age") +
  guides(shape = guide_legend(order = 1), col = guide_legend(order = 2))
ggsave(file="Graphs/Final_Socsim_Exp4_ASFR_ASMR.jpeg", width=17, height=9, dpi=200)


#----------------------------------------------------------------------------------------------------
## Summary measures: TFR and e0 ----
# Here, we use the rates by 1 year age group and 1 calendar year

# Total Fertility Rate ----
# Estimate Total Fertility Rate from asfr 1x1

# Estimate age-specific fertility rates 1x1 for the subset without 5% childless women
asfr_less_women_05_1 <- map_dfr(less_women_05, ~ estimate_fertility_rates(opop = .x,
                                                                          final_sim_year = 2022,
                                                                          year_min = 1750,
                                                                          year_max = 2023,
                                                                          year_group = 1,
                                                                          age_min_fert = 10,
                                                                          age_max_fert = 55, 
                                                                          age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_women_05_1, file = "Measures/asfr_less_women_05_1.RData")

# Estimate age-specific fertility rates 1x1 for the subset without 10% childless women
asfr_less_women_10_1 <- map_dfr(less_women_10, ~ estimate_fertility_rates(opop = .x,
                                                                          final_sim_year = 2022,
                                                                          year_min = 1750, 
                                                                          year_max = 2023, 
                                                                          year_group = 1,
                                                                          age_min_fert = 10, 
                                                                          age_max_fert = 55, 
                                                                          age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_women_10_1, file = "Measures/asfr_less_women_10_1.RData")

# Estimate age-specific fertility rates 1x1 for the subset without 20% childless women
asfr_less_women_20_1 <- map_dfr(less_women_20, ~ estimate_fertility_rates(opop = .x,
                                                                          final_sim_year = 2022, 
                                                                          year_min = 1750, 
                                                                          year_max = 2023, 
                                                                          year_group = 1,
                                                                          age_min_fert = 10, 
                                                                          age_max_fert = 55, 
                                                                          age_group = 1), 
                                     .id = "Sim_id") 
save(asfr_less_women_20_1, file = "Measures/asfr_less_women_20_1.RData")

# Load mean age-specific fertility rates 1x1 for the 10 simulations
load("Measures/asfr_whole_1.RData")
# Load asfr 1x1 for the subset with all direct ancestors and collateral kin
load("Measures/asfr_anc_zaukgausc_1.RData")
# Load asfr 1x1 for the subset without 5% childless women
load("Measures/asfr_less_women_05_1.RData")
# Load asfr 1x1 for the subset without 10% childless women
load("Measures/asfr_less_women_10_1.RData")
# Load asfr 1x1 for the subset without 20% childless women
load("Measures/asfr_less_women_20_1.RData")

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

# All Direct ancestors and collateral kin # Not included
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
  

# All Direct ancestors and collateral kin without 5% childless women
TFR_less_women_05 <- asfr_less_women_05_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "05% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 10% childless women
TFR_less_women_10 <- asfr_less_women_10_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "10% Omission",
         Rate = "TFR",           
         sex = "female")

# All Direct ancestors and collateral kin without 20% childless women
TFR_less_women_20 <- asfr_less_women_20_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "20% Omission",
         Rate = "TFR",           
         sex = "female")


## Plot TFR from whole SOCSIM simulation and subsets with different proportions of omitted children

bind_rows(TFR_whole, TFR_anc_zaukgausc,
          TFR_less_women_05, TFR_less_women_10, TFR_less_women_20) %>% 
  filter(Year >= 1850) %>%
  ggplot(aes(x = Year, y = TFR, group = Dataset, colour = Dataset)) +
  geom_line(linewidth = 1.3)+ 
  scale_color_viridis(option = "D", discrete = T, direction = -1)+
  theme_graphs() 
# labs(title = "Total Fertility Rate in Sweden (1751-2022), 
#retrieved from a SOCSIM simulation and subsets with different proportions of omitted women") 
ggsave(file="Graphs/socsim_Exp4_TFR.jpeg", width=17, height=9, dpi=200)


# Life Expectancy at birth ----
# Estimate life expectancy at birth from asmr 1x1 for the genealogical subsets 
# with different proportions of omitted childless women

# Estimate age-specific mortality rates for the subset without 5% childless women
asmr_less_women_05_1 <- map_dfr(less_women_05, ~ estimate_mortality_rates(opop = .x,
                                                                          final_sim_year = 2022, 
                                                                          year_min = 1750, 
                                                                          year_max = 2023, 
                                                                          year_group = 1,
                                                                          age_max_mort = 110, 
                                                                          age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_women_05_1, file = "Measures/asmr_less_women_05_1.RData")

# Calculate the mean and compute the life table
asmr_less_women_05_1 <- asmr_less_women_05_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
lt_less_women_05 <- lt_socsim(asmr_less_women_05_1)
save(lt_less_women_05, file = "Measures/lt_less_women_05.RData")


# Estimate age-specific mortality rates for the subset without 10% childless women
asmr_less_women_10_1 <- map_dfr(less_women_10, ~ estimate_mortality_rates(opop = .x,
                                                                          final_sim_year = 2022,
                                                                          year_min = 1750, 
                                                                          year_max = 2023, 
                                                                          year_group = 1,
                                                                          age_max_mort = 110,
                                                                          age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_women_10_1, file = "Measures/asmr_less_women_10_1.RData")

# Calculate the mean and compute the life table
asmr_less_women_10_1 <- asmr_less_women_10_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
lt_less_women_10 <- lt_socsim(asmr_less_women_10_1)
save(lt_less_women_10, file = "Measures/lt_less_women_10.RData")


# Estimate age-specific mortality rates for the subset without 20% childless women
asmr_less_women_20_1 <- map_dfr(less_women_20, ~ estimate_mortality_rates(opop = .x,
                                                                          final_sim_year = 2022, 
                                                                          year_min = 1750,
                                                                          year_max = 2023, 
                                                                          year_group = 1,
                                                                          age_max_mort = 110, 
                                                                          age_group = 1),
                                     .id = "Sim_id") 
save(asmr_less_women_20_1, file = "Measures/asmr_less_women_20_1.RData")

# Calculate the mean and compute the life table
asmr_less_women_20_1 <- asmr_less_women_20_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
lt_less_women_20 <- lt_socsim(asmr_less_women_20_1)
save(lt_less_women_20, file = "Measures/lt_less_women_20.RData")

# Wrangle data for plotting

# Load mean age-specific mortality rates 1x1 for the 10 simulations
load("Measures/asmr_whole_1.RData")
# Load life tables for the mean asmr for the 10 simulations
load("Measures/lt_whole.RData")
# Load life tables for the mean asmr from subset of all direct ancestors and collateral kin
load("Measures/lt_anc_zaukgausc.RData")
# Load life tables for the mean asmr from subset without 5% childless women
load("Measures/lt_less_women_05.RData")
# Load life tables for the mean asmr from subset without 10% childless women
load("Measures/lt_less_women_10.RData")
# Load life tables for the mean asmr from subset without 20% childless women
load("Measures/lt_less_women_20.RData")

# Year breaks. Extract all the unique numbers from the intervals 
year_breaks_mort_1 <- unique(as.numeric(str_extract_all(asmr_whole_1$year, "\\d+", simplify = T)))

# Year range to filter data . Check if necessary
year_range_mort_1 <- min(year_breaks_mort_1):max(year_breaks_mort_1-1)

# Whole SOCSIM simulation
lt_whole2 <- lt_whole %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Whole Simulation", 
         Rate = "e0") %>%    
  select(Year, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin 
lt_anc_zaukgausc2 <- lt_anc_zaukgausc %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct Ancestors + Collateral Kin",
         Rate = "e0") %>%    
  select(Year, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 5% childless women
lt_less_women_05b <- lt_less_women_05 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "05% Omission",
         Rate = "e0") %>%
  select(Year, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 10% childless women
lt_less_women_10b <- lt_less_women_10 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "10% Omission",
         Rate = "e0") %>%    
  select(Year, ex, Dataset, Rate, sex, Age)

# All Direct ancestors and collateral kin without 20% childless women
lt_less_women_20b <- lt_less_women_20 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "20% Omission",
         Rate = "e0") %>%    
  select(Year, ex, Dataset, Rate, sex, Age)

## Plot
bind_rows(lt_whole2, lt_anc_zaukgausc2,
          lt_less_women_05b, lt_less_women_10b, lt_less_women_20b) %>% 
  filter(Year >= 1850 & Age == 0) %>%
  ggplot(aes(x = Year, y = ex, group = Dataset, colour = Dataset)) +
  facet_grid(. ~ sex) +
  geom_line(linewidth = 1.3)+ 
  scale_color_viridis(option = "D", discrete = T, direction = -1)+
  theme_graphs() 
# labs(title = "Total Fertility Rate in Sweden (1751-2022), 
#retrieved from a SOCSIM simulation and subsets with different proportions of omitted women") 
ggsave(file="Graphs/socsim_Exp4_e0.jpeg", width=17, height=9, dpi=200)


#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

yrs_plot2 <- c(1850, 1900, 1950, 2000)

## TFR and e0 (for females) from whole SOCSIM simulation and subsets of "direct" and all collateral kin
bind_rows(TFR_whole %>% rename(Estimate = TFR),
          TFR_anc_zaukgausc %>% rename(Estimate = TFR),
          TFR_less_women_05 %>% rename(Estimate = TFR),
          TFR_less_women_10 %>% rename(Estimate = TFR),
          TFR_less_women_20 %>% rename(Estimate = TFR)) %>%
  mutate(sex = "female") %>%  
  bind_rows(lt_whole2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_anc_zaukgausc2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_less_women_05b %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_less_women_10b %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_less_women_20b %>% rename(Estimate = ex) %>% filter(Age == 0)) %>%
  filter(sex == "female" & Year >= 1850) %>%
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life Expectancy at Birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life Expectancy at Birth"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), 
             aes(shape = Dataset), size = 11) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#FF5D00", "#8C463B", "#8700AA", "#1A3077", "#1A7761"))+
  scale_shape_manual(values = c(21, 22, 23, 18, 46)) +
  #scale_x_continuous(breaks = yrs_plot)+
  theme_graphs()
# labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and subsets with different proportions of omitted children")
# Save the plot
ggsave(file="Graphs/Final_Socsim_Exp4_TFR_e0.jpeg", width=17, height=9, dpi=200)