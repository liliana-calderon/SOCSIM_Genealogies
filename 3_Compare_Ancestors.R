#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Genealogies - Trace and compare subset of direct ancestors of a given ego(s) 
# U:/SOCSIM/SOCSIM_Genealogies/3_Compare_Ancestors.R

## Trace direct ancestors of a given ego(s) from a SOCSIM microsimulation for Sweden (1751-2022)
# and compare demographic measures from the whole simulation and the genealogical subsets

# Created by Liliana Calderon on 23-09-2022
# Last modified by Liliana Calderon on 31-05-2023

## NB: To run this code, it is necessary to have already run the script 1_Run_Simulations.R

#------------------------------------------------------------------------------------------------------
## General settings ----

# Prevent scientific notation (useful for the rate calculation)
options(scipen=999999)

## Load packages 
library(tidyverse)
library(ggh4x)  # To facet scales-
library(svglite) # To save svg files
library(viridis)
library(rsocsim) # Functions to estimate rates

## Load functions to get direct ancestors
# source("Functions_Ancestors.R")
source("Functions_Ancestors_Imp.R") # Improved function

## Load functions to estimate age-specific fertility rates 
# This is a slightly modified version of the functions in the package 
# that allows to handle the intentional duplicates in the data
source("Functions_Fertility_Rates_Mod.R")

## Load theme for the graphs
source("Functions_Graphs.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with asmr calculated with the estimate_mortality_rates
source("Functions_Life_Table.R")

#------------------------------------------------------------------------------------------------------
## Read the output .opop file ----

## Randomly choose the simulation seed to use 
load("sims_seeds.rda")
# seed <-  sample(sims_seeds, 1, replace = F) 
seed <- "1129"

## We use only one of the 10 simulations, same seed chosen in 3_Compare ancestors
opop <- read_opop(folder = getwd(), supfile = "Sweden.sup", seed = seed, 
                  suffix = "",  fn = NULL)

#------------------------------------------------------------------------------------------------------
## Trace direct ancestors of people alive in 2022 as a proxy of current genealogists 

# Pids of people alive at the end of the simulation, i.e. dod == 0, 
# who are older than 18 years old on 01-01-2023, i,e. dob 1914-2004
egos2023 <- opop %>% 
  mutate(last_month = max(dob),
         final_sim_year = 2022, ## Change if necessary
         Generation = asYr(dob, last_month, final_sim_year)) %>%  
  filter(dod == 0 & Generation <= final_sim_year-18) %>% 
  pull(pid)

# Get a sample of 10% of people alive in 2022. 
sample_size <- round(length(egos2023)/10)
egos2023_samp <-  sample(egos2023, sample_size, replace = F)
save(egos2023_samp, file = "egos2023_samp_10.RData")

## Map the function to get the ancestors of a sample of individuals alive in 2022 (older than 18 years)
start <- Sys.time()
ancestors_egos2023_10 <- map_dfr(egos2023_samp, get_ancestors) %>%
  left_join(select(opop, c(pid, fem, dob, dod, mom, marid, mstat)), by = "pid")
end <- Sys.time()
print(end-start)
# Time difference of 10.57319 hours for 44094 egos
# Before nearly 70 hours for 17832 egos

# Save the data frame
save(ancestors_egos2023_10, file = "ancestors_egos2023_10.RData")

#----------------------------------------------------------------------------------------------------
## Recover age-specific fertility and mortality rates  -----
# Retrieve and compare the rates derived from the whole simulation with those from the genealogical subset
# of direct ancestors both will and without duplicates

## Calculate ASFR and ASMR for the Whole simulation (seed "13486")

# Retrieve age-specific fertility rates for the whole single simulation 
asfr_whole <- estimate_fertility_rates(opop = opop,
                                      final_sim_year = 2022 , #[Jan-Dec]
                                      year_min = 1750, # Closed [
                                      year_max = 2020, # Open )
                                      year_group = 5, 
                                      age_min_fert = 10, # Closed [
                                      age_max_fert = 55, # Open )
                                      age_group = 5) #[,)
save(asfr_whole, file = "Measures/asfr_whole.RData")

# Retrieve age-specific mortality rates for the whole single simulation
asmr_whole <- estimate_mortality_rates(opop = opop, 
                                      final_sim_year = 2022, #[Jan-Dec]
                                      year_min = 1750, # Closed
                                      year_max = 2020, # Open )
                                      year_group = 5,
                                      age_max_mort = 110, # Open )
                                      age_group = 5) #[,)
save(asmr_whole, file = "Measures/asmr_whole.RData")

# Load the data frame with the ancestors of 10% sample of egos alive in 2022
load("ancestors_egos2023_10.RData")

#  Calculate ASFR and ASMR for genealogical subset of direct ancestors of population alive in 01-01-2022 with duplicates 

# Copy the vector of direct ancestors with duplicates. 
ancestors_egos2023_wd <- ancestors_egos2023_10 

# Retrieve age-specific fertility rates for the genealogical subset of direct ancestors with duplicates
asfr_wd <- estimate_fertility_rates_mod(opop = ancestors_egos2023_wd,
                                        final_sim_year = 2022 , #[Jan-Dec]
                                        year_min = 1750, # Closed [
                                        year_max = 2020, # Open )
                                        year_group = 5, 
                                        age_min_fert = 10, # Closed [
                                        age_max_fert = 55, # Open )
                                        age_group = 5) #[,)
save(asfr_wd, file = "Measures/asfr_wd.RData")

# Retrieve age-specific mortality rates for the genealogical subset of direct ancestor with duplicates
asmr_wd <- estimate_mortality_rates(opop = ancestors_egos2023_wd,
                                   final_sim_year = 2022, #[Jan-Dec]
                                   year_min = 1750, # Closed
                                   year_max = 2020, # Open )
                                   year_group = 5,
                                   age_max_mort = 110, # Open )
                                   age_group = 5) #[,)
save(asmr_wd, file = "Measures/asmr_wd.RData")


## Calculate ASFR and ASMR for genealogical subset of direct ancestors without duplicates 

# Ancestors without duplicates for sample 
ancestors_egos2023_wod <- ancestors_egos2023_10 %>% distinct(pid, .keep_all = TRUE)

# Retrieve age-specific fertility rates for the genealogical subset of direct ancestors without duplicates
asfr_wod <- estimate_fertility_rates(opop = ancestors_egos2023_wod,
                                    final_sim_year = 2022 , #[Jan-Dec]
                                    year_min = 1750, # Closed [
                                    year_max = 2020, # Open )
                                    year_group = 5, 
                                    age_min_fert = 10, # Closed [
                                    age_max_fert = 55, # Open )
                                    age_group = 5) #[,)
save(asfr_wod, file = "Measures/asfr_wod.RData")

# Retrieve age-specific mortality rates for the genealogical subset of direct ancestors without duplicates
asmr_wod <- estimate_mortality_rates(opop = ancestors_egos2023_wod,
                                      final_sim_year = 2022, #[Jan-Dec]
                                      year_min = 1750, # Closed
                                      year_max = 2020, # Open )
                                      year_group = 5,
                                      age_max_mort = 110, # Open )
                                      age_group = 5) #[,)
save(asmr_wod, file = "Measures/asmr_wod.RData")

#----------------------------------------------------------------------------------------------------
## Plot results for the genealogical subsets of direct ancestors with and without duplicates ----

# Load ASFR and ASMR for the genealogical subset of direct ancestors with duplicates
load("Measures/asfr_wd.RData")
load("Measures/asmr_wd.RData")

# Load ASFR and ASMR for the genealogical subset of direct ancestors without duplicates
load("Measures/asfr_wod.RData")
load("Measures/asmr_wod.RData")

# Choose years to plot (in intervals).
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_wd$age)

## ASFR and ASMR (for women) genealogical subset of direct ancestors with duplicates
bind_rows(asfr_wd %>% 
            mutate(rate = "ASFR",
                   sex = "female"),
          asmr_wd %>% 
            mutate(rate = "ASMR") %>% 
            filter(sex == "female")) %>% 
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(socsim !=0 & !is.infinite(socsim) & !is.nan(socsim)) %>% 
  filter(year %in% yrs_plot) %>% 
  ggplot(aes(x = age, y = socsim, group = year, colour = year)) +
  geom_line(linewidth = 1) +
  facet_wrap(. ~ rate, scales = "free") + 
  theme_graphs()  +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10"))) +
  scale_color_viridis(option = "F", discrete = T, direction = -1)+
  labs(x = "Age", y = "Estimate")


## ASFR and ASMR (for women) genealogical subset of direct ancestors without duplicates
bind_rows(asfr_wod %>% 
            mutate(rate = "ASFR",
                   sex = "female"),
          asmr_wod %>% 
            mutate(rate = "ASMR") %>% 
            filter(sex == "female")) %>% 
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(socsim !=0 & !is.infinite(socsim) & !is.nan(socsim)) %>% 
  filter(year %in% yrs_plot) %>% 
  ggplot(aes(x = age, y = socsim, group = year, colour = year)) +
  geom_line(linewidth = 1) +
  facet_wrap(. ~ rate, scales = "free") + 
  theme_graphs()  +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10"))) +
  scale_color_viridis(option = "F", discrete = T, direction = -1)+
  labs(x = "Age", y = "Estimate")


#----------------------------------------------------------------------------------------------------
## Comparison of a whole SOCSIM simulation with some genealogical subsets of direct ancestors ----

# Load ASFR and ASMR for the whole single simulation 
load("Measures/asfr_whole.RData")
load("Measures/asmr_whole.RData")

#### Age-specific Fertility Rates ----

# Whole SOCSIM simulation
asfr_whole2 <- asfr_whole %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Whole_simulation", 
         Rate = "ASFR") 

# Genealogical subset of direct ancestors with duplicates
asfr_wd2 <- asfr_wd %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Ancestors_w_dup",
         Rate = "ASFR") 

# Genealogical subset of direct ancestors without duplicates
asfr_wod2 <- asfr_wod %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Ancestors_wo_dup", 
         Rate = "ASFR") 

## Plot ASFR from whole SOCSIM simulation and the genealogical subset of direct ancestors

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2, asfr_wd2, asfr_wod2) %>% 
  filter(year %in% yrs_plot) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset)))+
  geom_line(aes(colour = year, linetype = Dataset), linewidth = 1.2)+ 
  scale_color_viridis(option = "D", discrete = T, direction = -1) +
  scale_linetype_manual(values = c("11", "22", "solid")) +
  theme_graphs()
# labs(title = "Age-Specific Fertility Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and subsets of direct ancestors") 
ggsave(file="Graphs/Socsim_Exp1_ASFR.jpeg", width=17, height=9, dpi=200)


## Age-Specific Mortality rates ----

# Whole SOCSIM simulation
asmr_whole2 <- asmr_whole %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Whole_simulation",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)
  
# Genealogical subset of direct ancestors with duplicates
asmr_wd2 <- asmr_wd %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Ancestors_w_dup",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Genealogical subset of direct ancestors without duplicates
asmr_wod2 <- asmr_wod %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Ancestors_wo_dup",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)


## Plotting ASMR from whole SOCSIM simulation and some genealogical subsets

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asmr_whole2, asmr_wd2, asmr_wod2) %>% 
  rename(Year = year) %>% 
  filter(Year %in% yrs_plot) %>%  
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>% 
  ggplot(aes(x = age, y = mx, group = interaction(Year, Dataset)))+
  facet_wrap(~Sex) +
  geom_line(aes(colour = Year, linetype = Dataset), linewidth = 1.2)+ 
  scale_color_viridis(option = "G", discrete = T, direction = -1) +
  scale_linetype_manual(values = c("11", "22", "solid")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs()
  #labs(title = "Age-Specific Mortality Rates in Sweden (1751-2022), retrieved from a SOCSIM simulation and subset of direct ancestors") 
ggsave(file="Graphs/socsim_Exp1_ASMR.jpeg", width=17, height=9, dpi=200)


## Final plot combining ASFR and ASMR ----

# Change years to plot only to two periods
yrs_plot <- c("[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_whole2$age)

## Ploting ASFR and ASMR (for females) from whole SOCSIM simulation and some genealogical subsets
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_wd2 %>% rename(Estimate = ASFR),
          asfr_wod2 %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx),
            asmr_wd2 %>% rename(Estimate = mx),
            asmr_wod2 %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels), 
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates"),
         Dataset = case_when(Dataset == "Ancestors_w_dup" ~ "Experiment 1 with duplicates",
                             Dataset == "Ancestors_wo_dup" ~ "Experiment 1 without duplicates",
                             TRUE ~ "Whole simulation")) %>% 
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(linewidth = 1.3)+
  geom_point(aes(shape = Dataset), size = 7)+
  scale_color_manual(values = c("#B72779", "#2779B7"))+
  scale_shape_manual(values = c(15,17,46)) + 
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(x = "Age")
ggsave(file="Graphs/Final_Socsim_Exp1_ASFR_ASMR.jpeg", width=17, height=9, dpi=200)

## Save as .svg file for poster
# ggsave(file="Graphs/Socsim_Exp1_ASFR_ASMR.svg", device = "svg", units = "in", width=15, height=8, dpi=200) 

#----------------------------------------------------------------------------------------------------
#### Summary measures: TFR and e0 ----
# Here, we use the rates by single calendar year and 1 year of age

#### Total Fertility Rate ----
# Calculate Total Fertility Rate from asfr 1x1

# Retrieve age-specific fertility rates 1x1 for the whole single simulation 
asfr_whole_1 <- estimate_fertility_rates(opop = opop,
                                final_sim_year = 2022 , #[Jan-Dec]
                                year_min = 1750, # Closed [
                                year_max = 2020, # Open )
                                year_group = 1, 
                                age_min_fert = 10, # Closed [
                                age_max_fert = 55, # Open )
                                age_group = 1) #[,)
save(asfr_whole_1, file = "Measures/asfr_whole_1.RData")

# Retrieve age-specific fertility rates 1x1 for the genealogical subset of direct ancestors with duplicates
asfr_wd_1 <- estimate_fertility_rates_mod(opop = ancestors_egos2023_wd,
                                          final_sim_year = 2022 , #[Jan-Dec]
                                          year_min = 1750, # Closed [
                                          year_max = 2020, # Open )
                                          year_group = 1, 
                                          age_min_fert = 10, # Closed [
                                          age_max_fert = 55, # Open )
                                          age_group = 1) #[,)
save(asfr_wd_1, file = "Measures/asfr_wd_1.RData")

# Retrieve age-specific fertility rates 1x1 for the genealogical subset of direct ancestors without duplicates
asfr_wod_1 <- estimate_fertility_rates(opop = ancestors_egos2023_wod,
                                      final_sim_year = 2022 , #[Jan-Dec]
                                      year_min = 1750, # Closed [
                                      year_max = 2020, # Open )
                                      year_group = 1, 
                                      age_min_fert = 10, # Closed [
                                      age_max_fert = 55, # Open )
                                      age_group = 1) #[,)
save(asfr_wod_1, file = "Measures/asfr_wod_1.RData")

## Load asfr_1 data
load("Measures/asfr_whole_1.RData")
load("Measures/asfr_wd_1.RData")
load("Measures/asfr_wod_1.RData")

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
  mutate(Dataset = "Whole_simulation",
         Rate = "TFR", 
         sex = "female")

# Genealogical subset of direct ancestors with duplicates
TFR_wd <- asfr_wd_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "Ancestors_w_dup",
         Rate = "TFR", 
         sex = "female") 

# Genealogical subset of direct ancestors without duplicates
TFR_wod <- asfr_wod_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "Ancestors_wo_dup",
         Rate = "TFR", 
         sex = "female")

## Plott TFR from whole SOCSIM simulation and genealogical subsets of direct ancestors with(out) duplicates
bind_rows(TFR_whole, TFR_wd, TFR_wod) %>% 
  filter(Year >= 1751) %>%
  ggplot(aes(x = Year, y = TFR)) +
  geom_line(aes(colour = Dataset, linetype = Dataset), linewidth = 1.3)+ 
  scale_color_manual(values = c("#FC8961", "#B73779", "#1A7761"))+
  scale_linetype_manual(values = c("11", "22", "solid")) +
  theme_graphs() 
# labs(title = "Total Fertility Rate in Sweden (1751-2022), retrieved from a SOCSIM simulation and genealogical subset of direct ancestors") 
ggsave(file="Graphs/socsim_Exp1_TFR.jpeg", width=17, height=9, dpi=200)

## Life Expectancy at birth ----
# Calculate life expectancy at birth from asmr 1x1

# Retrieve age-specific mortality rates for the whole simulation
asmr_whole_1 <- estimate_mortality_rates(opop = opop, 
                                        final_sim_year = 2022, #[Jan-Dec]
                                        year_min = 1750, # Closed
                                        year_max = 2020, # Open )
                                        year_group = 1,
                                        age_max_mort = 110, # Open )
                                        age_group = 1) #[,)
save(asmr_whole_1, file = "Measures/asmr_whole_1.RData")


# Retrieve age-specific mortality rates for the genealogical subset of direct ancestors with duplicates
asmr_wd_1 <- estimate_mortality_rates(opop = ancestors_egos2023_wd,
                                     final_sim_year = 2022, #[Jan-Dec]
                                     year_min = 1750, # Closed
                                     year_max = 2020, # Open )
                                     year_group = 1,
                                     age_max_mort = 110, # Open )
                                     age_group = 1) #[,)
save(asmr_wd_1, file = "Measures/asmr_wd_1.RData")

# Retrieve age-specific mortality rates for the genealogical subset of direct ancestors without duplicates
asmr_wod_1 <- estimate_mortality_rates(opop = ancestors_egos2023_wod,
                                      final_sim_year = 2022, #[Jan-Dec]
                                      year_min = 1750, # Closed
                                      year_max = 2020, # Open )
                                      year_group = 1,
                                      age_max_mort = 110, # Open )
                                      age_group = 1) #[,)
save(asmr_wod_1, file = "Measures/asmr_wod_1.RData")

# Compute life table from asmr 1x1 for Whole SOCSIM simulation
lt_whole <- lt_socsim(asmr_whole_1)
save(lt_whole, file = "Measures/lt_whole.RData")

# Compute life table from asmr 1x1 for the genealogical subset of direct ancestors with duplicates
lt_wd <- lt_socsim(asmr_wd_1)
save(lt_wd, file = "Measures/lt_wd.RData")

# Compute life table from asmr 1x1 for the genealogical subset of direct ancestors without duplicates
lt_wod <- lt_socsim(asmr_wod_1)
save(lt_wod, file = "Measures/lt_wod.RData")


## Load lt_1 data
load("Measures/asmr_whole_1.RData")
load("Measures/lt_whole.RData")
load("Measures/lt_wd.RData")
load("Measures/lt_wod.RData")

# Year breaks. Extract all the unique numbers from the intervals 
year_breaks_mort_1 <- unique(as.numeric(str_extract_all(asmr_whole_1$year, "\\d+", simplify = T)))

# Year range to filter HMD data
year_range_mort_1 <- min(year_breaks_mort_1):max(year_breaks_mort_1-1)

# Whole SOCSIM simulation
lt_whole2 <- lt_whole %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Whole_simulation",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Genealogical subset of direct ancestors with duplicates
lt_wd2 <- lt_wd %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Ancestors_w_dup",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Genealogical subset of direct ancestors without duplicates
lt_wod2 <- lt_wod %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Ancestors_wo_dup",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

bind_rows(lt_whole2, lt_wd2, lt_wod2) %>% 
  filter(Age == 0 & Year %in% year_range_mort_1) %>%
  mutate(Sex = ifelse(sex == "female", "Female", "Male")) %>% 
  ggplot(aes(x = Year, y = ex, group = Dataset))+
  geom_line(aes(colour = Dataset, linetype = Dataset), linewidth = 1.3)+
  scale_color_manual(values = c("#FC8961", "#B73779", "#51127C"))+
  scale_linetype_manual(values = c("11", "22", "solid")) +
  facet_wrap(~Sex) +
  theme_graphs()+
  labs(title = "Life expectancy at birth in Sweden (e0), 1751-2020, a SOCSIM simulation and genealogical subset of direct ancestors",
       y = "e0") 
ggsave(file="Graphs/socsim_Exp1_e0.jpeg", width=17, height=9, dpi=200)


#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

## TFR and e0 (for females) from whole SOCSIM simulation and genealogical subsets of direct ancestors

yrs_plot2 <- c(1750, 1800, 1850, 1900, 1950, 2000)

bind_rows(TFR_whole %>% 
          rename(Estimate = TFR), 
          TFR_wd %>% 
          rename(Estimate = TFR), 
          TFR_wod %>%           
          rename(Estimate = TFR)) %>%  
  bind_rows(lt_whole2 %>% 
              rename(Estimate = ex) %>% 
              filter(Age == 0),
            lt_wd2 %>% 
              rename(Estimate = ex) %>% 
              filter(Age == 0),
            lt_wod2 %>% 
              rename(Estimate = ex) %>% 
              filter(Age == 0)) %>% 
  filter(sex == "female") %>% 
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life expectancy at birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life expectancy at birth")),
         Dataset = case_when(Dataset == "Ancestors_w_dup" ~ "Experiment 1 with duplicates",
                             Dataset == "Ancestors_wo_dup" ~ "Experiment 1 without duplicates",
                             TRUE ~ "Whole simulation")) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), aes(shape = Dataset), size = 9)+
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#771A30", "#331A77", "#1A7761"))+
  scale_shape_manual(values = c(15,17,46)) + 
  scale_x_continuous(breaks = yrs_plot2)+
  theme_graphs()
# labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2020), retrieved from HFD, HMD and 10 SOCSIM simulation outputs") + 

# Save the plot
ggsave(file="Graphs/Final_Socsim_Exp1_TFR_e0.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Sex Ratio at Birth and Infant Mortality Rate
# The Functions_Retrieve_Rates.R must be called to use the asYr() function

## Define years of not set in the Global Environment
final_sim_year <- 2022 #[Jan-Dec]
year_min <- 1750 # Closed [
year_max <- 2020 # Open )

# Year range
year_range <- year_min:(year_max-1)

# Find last month of the simulation
last_month <- max(opop$dob)

# Sex Ratio at Birth by year for the whole simulation
SRB_whole <- opop %>% 
  mutate(Year = asYr(dob, last_month, final_sim_year),
         Sex = ifelse(fem == 1, "Female", "Male")) %>% 
  filter(Year %in% year_range) %>% 
  count(Year, Sex) %>%
  mutate(Year = factor(Year, levels = year_range), 
         Sex = factor(Sex, levels = c("Female", "Male"))) %>% 
  complete(Year, Sex, fill = list(n = 0)) %>% 
  mutate(Year = as.numeric(as.character(Year)),
         Sex = as.character(Sex)) %>% 
  pivot_wider(names_from = Sex, values_from = n) %>% 
  mutate(SRB = Male/Female, 
         Measure = "SRB",
         Dataset = "Whole_simulation") 

# Sex Ratio at Birth by year for genealogical subset of direct ancestors with duplicates
SRB_wd <- ancestors_egos2023_wd %>% 
  mutate(Year = asYr(dob, last_month, final_sim_year),
         Sex = ifelse(fem == 1, "Female", "Male")) %>% 
  filter(Year %in% year_range) %>% 
  count(Year, Sex) %>%
  mutate(Year = factor(Year, levels = year_range), 
         Sex = factor(Sex, levels = c("Female", "Male"))) %>% 
  complete(Year, Sex, fill = list(n = 0)) %>% 
  mutate(Year = as.numeric(as.character(Year)),
         Sex = as.character(Sex)) %>% 
  pivot_wider(names_from = Sex, values_from = n) %>% 
  mutate(SRB = Male/Female,         
         Measure = "SRB",
         Dataset = "Ancestors_w_dup") 

# Sex Ratio at Birth by year for genealogical subset of direct ancestors without duplicates
SRB_wod <- ancestors_egos2023_wod %>% 
  mutate(Year = asYr(dob, last_month, final_sim_year),
         Sex = ifelse(fem == 1, "Female", "Male")) %>% 
  filter(Year %in% year_range) %>% 
  count(Year, Sex) %>%
  mutate(Year = factor(Year, levels = year_range), 
         Sex = factor(Sex, levels = c("Female", "Male"))) %>% 
  complete(Year, Sex, fill = list(n = 0)) %>% 
  mutate(Year = as.numeric(as.character(Year)),
         Sex = as.character(Sex)) %>% 
  pivot_wider(names_from = Sex, values_from = n) %>% 
  mutate(SRB = Male/Female,         
         Measure = "SRB",
         Dataset = "Ancestors_wo_dup") 
  
# Plotting SRB
bind_rows(SRB_whole, SRB_wd, SRB_wod) %>% 
filter(Year >= 1751 & !is.na(SRB)) %>% # No births after 2003
  ggplot(aes(x = Year, y = SRB, group = Dataset, color = Dataset, shape = Dataset))+
  geom_point(data = . %>% filter(Year %in% yrs_plot2), size = 9)+
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#771A30", "#331A77", "#1A7761"))+
  scale_shape_manual(values = c(15,17,46)) + 
  scale_x_continuous(breaks = yrs_plot2)+
  theme_graphs()

#### Infant Mortality Rate, both sexes

# Births by year from the whole simulation
Births_whole <- opop %>% 
    mutate(Year = asYr(dob, last_month, final_sim_year)) %>% 
    filter(Year %in% year_range) %>% 
    count(Year) %>%
    mutate(Year = factor(Year, levels = year_range)) %>% 
    complete(Year, fill = list(n = 0))  %>%
    mutate(Year = as.numeric(as.character(Year)),
           Dataset = "Whole_simulation", 
           Event = "Births")
  
# Births by year from genealogical subset of direct ancestors with duplicates
Births_wd <- ancestors_egos2023_wd %>% 
    mutate(Year = asYr(dob, last_month, final_sim_year)) %>% 
    filter(Year %in% year_range) %>% 
    count(Year) %>%
    mutate(Year = factor(Year, levels = year_range)) %>% 
    complete(Year, fill = list(n = 0))  %>%
    mutate(Year = as.numeric(as.character(Year)),
           Dataset = "Ancestors_w_dup", 
           Event = "Births")
  
# Births by year from genealogical subset of direct ancestors without duplicates
Births_wod <- ancestors_egos2023_wod %>% 
    mutate(Year = asYr(dob, last_month, final_sim_year)) %>% 
    filter(Year %in% year_range) %>% 
    count(Year) %>%
    mutate(Year = factor(Year, levels = year_range)) %>% 
    complete(Year, fill = list(n = 0))  %>%
    mutate(Year = as.numeric(as.character(Year)),
           Dataset = "Ancestors_wo_dup", 
           Event = "Births")
  
# Deaths below age 1 (0-11 months) by year from the whole simulation
Deaths_0_whole <- opop %>% 
  filter(dod != 0) %>% 
  mutate(age_death_months = dod-dob,
         Year = asYr(dod, last_month, final_sim_year)) %>% 
  filter(age_death_months < 12 & Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>% 
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Whole_simulation", 
         Event = "Deaths")

# Deaths below age 1 (0-11 months) from genealogical subset of direct ancestors with duplicates
# There should be no infant mortality in this subset, but let's double check it
Deaths_0_wd <- ancestors_egos2023_wd %>% 
  filter(dod != 0) %>% 
  mutate(age_death_months = dod-dob,
         Year = asYr(dod, last_month, final_sim_year)) %>% 
  filter(age_death_months < 12 & Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>% 
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Ancestors_w_dup", 
         Event = "Deaths")

# Deaths below age 1 (0-11 months) from genealogical subset of direct ancestors without duplicates
# There should be no infant mortality in this subset, but let's double check it
Deaths_0_wod <- ancestors_egos2023_wod %>% 
  filter(dod != 0) %>% 
  mutate(age_death_months = dod-dob,
         Year = asYr(dod, last_month, final_sim_year)) %>% 
  filter(age_death_months < 12 & Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>% 
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Ancestors_wo_dup", 
         Event = "Deaths")

# Calculate and Plot Infant Mortality Rate (IMR)
IMR <- bind_rows(Births_whole, Births_wd, Births_wod, Deaths_0_whole, Deaths_0_wd, Deaths_0_wod) %>%
  pivot_wider(names_from = Event, values_from = n) %>% 
  mutate(IMR = Deaths/Births,
         Measure = "IMR") 

IMR %>% 
  filter(Year >= 1751) %>% 
  ggplot(aes(x = Year, y = IMR, group = Dataset, color = Dataset))+
  geom_line(linewidth = 1)+
  scale_color_manual(values = c("#771A30", "#331A77", "#1A7761"))+
  theme_graphs() 
# There is no IMR for the genealogical subsets

# Plot SRB and IMR together
bind_rows(SRB_whole, SRB_wd, SRB_wod) %>% 
  select(Year, Dataset, SRB) %>% 
  full_join(IMR %>% select(Year, Dataset, IMR), by = c("Year", "Dataset")) %>% 
  pivot_longer(SRB:IMR, names_to = "Measure", values_to = "Value") %>% 
  filter(Year >= 1751 & !is.na(Value)) %>%
  mutate(Measure = ifelse(Measure == "SRB", "Sex Ratio at Birth", "Infant Mortality Rate"), 
         Measure = factor(Measure, levels = c("Sex Ratio at Birth", "Infant Mortality Rate")),
         Dataset = case_when(Dataset == "Ancestors_w_dup" ~ "Experiment 1 with duplicates",
                             Dataset == "Ancestors_wo_dup" ~ "Experiment 1 without duplicates",
                             TRUE ~ "Whole simulation")) %>%
  ggplot(aes(x = Year, y = Value, group = Dataset, color = Dataset, shape = Dataset))+
  facet_wrap(. ~ Measure, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), aes(shape = Dataset), size = 9)+
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#771A30", "#331A77", "#1A7761"))+
  scale_shape_manual(values = c(15,17,46)) +
  theme_graphs() +
  facetted_pos_scales(y = list(SRB = scale_y_continuous(limits=c(0.7, 1.3)),
                               IMR =  scale_y_continuous())) +
  theme_graphs() +
  theme(axis.title.y = element_blank())
ggsave(file="Graphs/Socsim_Exp1_SRB_IMR.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Births and Deaths by year from whole simulation and genealogical subsets of direct ancestors -----

## final_sim_year, year_min, year_max, year_range must be set in the Global Environment. 
# They are defined above. 
# For the Birth counts we use the same df calculated before

# Death counts by year from the whole simulation
Deaths_whole <- opop %>% 
  filter(dod != 0) %>% 
  mutate(Year = asYr(dod, last_month, final_sim_year)) %>% 
  filter(Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>%
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Whole_simulation", 
         Event = "Deaths")

# Death counts by year from genealogical subset of direct ancestors with duplicates
Deaths_wd <- ancestors_egos2023_wd %>% 
  filter(dod != 0) %>% 
  mutate(Year = asYr(dod, last_month, final_sim_year)) %>% 
  filter(Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>%
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Ancestors_w_dup", 
         Event = "Deaths")

# Death counts by year from genealogical subset of direct ancestors without duplicates
Deaths_wod <- ancestors_egos2023_wod %>% 
  filter(dod != 0) %>% 
  mutate(Year = asYr(dod, last_month, final_sim_year)) %>% 
  filter(Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>%
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Ancestors_wo_dup", 
         Event = "Deaths")

# Plotting birth and death counts together. 
bind_rows(Births_whole, Births_wd, Births_wod, Deaths_whole, Deaths_wd, Deaths_wod) %>% 
  filter(Year >= 1751 & n!=0) %>%
  mutate(Dataset = case_when(Dataset == "Ancestors_w_dup" ~ "Experiment 1 with duplicates",
                             Dataset == "Ancestors_wo_dup" ~ "Experiment 1 without duplicates",
                             TRUE ~ "Whole simulation")) %>%
  ggplot(aes(x = Year, y = n, group = Dataset, color = Dataset, shape = Dataset))+
  facet_wrap(. ~ Event) + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), aes(shape = Dataset), size = 7)+
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#771A30", "#331A77", "#1A7761"))+ 
  scale_shape_manual(values = c(15,17,46)) +
  theme_graphs() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(file="Graphs/Socsim_Exp1_Births_Deaths.jpeg", width=17, height=9, dpi=200)