#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Genealogies - Trace and compare subset of direct ancestors of a given ego(s) 
# U:/SOCSIM/SOCSIM_Genealogies/3_Compare_Ancestors.R

## Trace direct ancestors of a given ego(s) from a SOCSIM microsimulation for Sweden (1751-2022)
# and compare demographic measures from the whole simulation and the genealogical subsets

# Created by Liliana Calderon on 23-09-2022
# Last modified by Liliana Calderon on 21-07-2023

## NB: To run this code, it is necessary to have already run the script 1_Run_Simulations.R
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

## Load function to get direct ancestors
source("Functions_Ancestors.R")

## Load function to estimate age-specific fertility rates 
# This is a slightly modified version of the function in the package 
# that allows to handle the intentional duplicates in the data
source("Functions_Fertility_Rates_Mod.R")

## Load theme for the graphs and to convert SOCSIM time
source("Functions_Graphs.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with asmr calculated with rsocsim::estimate_mortality_rates()
source("Functions_Life_Table.R")

#------------------------------------------------------------------------------------------------------
## Trace direct ancestors of people alive in 2023 ----

# Load saved list with opop from 10 simulations, generated in 1_Run_Simulations.R
load("sims_opop.RData")
# Load saved list with omar from 10 simulations, generated in 1_Run_Simulations.R
load("sims_omar.RData")
# Create a sub-folder called "Subsets" to save the opop subsets used for each experiment
ifelse(!dir.exists("Subsets"), dir.create("Subsets"), FALSE)

# Function to get a 10% sample of people alive at the end of each simulation, i.e. dod == 0, 
# who are older than 18 years old on 01-01-final_sim_year+1 

sample_egos <- function(opop = opop, final_sim_year, percentage) {
  egos <- opop %>% 
    mutate(last_month = max(dob),
           Generation = asYr(dob, last_month, final_sim_year)) %>%  
    filter(dod == 0 & Generation <= final_sim_year-18) %>% 
    pull(pid)
  sample_size <- round(length(egos)*percentage/100)
  egos_samp <-  sample(egos, sample_size, replace = F)
  return(egos_samp)
}

# Get sample for each of the 10 simulations
egos_samp_10 <- map(sims_opop, ~ sample_egos(opop = .x,
                                             final_sim_year = 2022, 
                                             percentage = 10)) 
# save(egos_samp_10, file = "Subsets/egos_samp_10.RData")

## Retrieve the ancestors of each simulation sample of egos alive in 2023
ancestors_10 <- map2_dfr(egos_samp_10, sims_opop,
                           ~ retrieve_ancestors(egos = .x, opop = .y), 
                           .id = "Sim_id") 
# Save the data frame with the ancestors of 10 simulations samples
save(ancestors_10, file = "Subsets/ancestors_10.RData")

#----------------------------------------------------------------------------------------------------
## Age-Specific Fertility and Mortality rates, 5x5  -----
# Estimate and compare the rates derived from the whole simulation with those from the genealogical subset
# of direct ancestors both with and without duplicates

# Load the data frame with the ancestors of the 10 simulations samples of egos alive in 2023
load("Subsets/ancestors_10.RData")

#  Estimate ASFR and ASMR for genealogical subset of direct ancestors of population alive in 01-01-2023 with duplicates 

# Create a list of data frames containing opop of ancestors for each simulation
ancestors_dir_wd <- ancestors_10 %>%
  split(.$Sim_id)

# Estimate age-specific fertility rates for the genealogical subset of direct ancestors with duplicates
# We need to use here the modified function that allows the joining of data frames with duplicates
asfr_dir_wd <- map_dfr(ancestors_dir_wd, ~ estimate_fertility_rates_mod(opop = .x,
                                                                        final_sim_year = 2022, #[Jan-Dec]
                                                                        year_min = 1750, # Closed [
                                                                        year_max = 2020, # Open )
                                                                        year_group = 5, 
                                                                        age_min_fert = 10, # Closed [
                                                                        age_max_fert = 55, # Open )
                                                                        age_group = 5), # [,)
                       .id = "Sim_id") 
save(asfr_dir_wd, file = "Measures/asfr_dir_wd.RData")

# Estimate age-specific mortality rates for the genealogical subset of direct ancestor with duplicates
asmr_dir_wd <- map_dfr(ancestors_dir_wd, ~ estimate_mortality_rates(opop = .x,
                                                                    final_sim_year = 2022, #[Jan-Dec]
                                                                    year_min = 1750, # Closed
                                                                    year_max = 2020, # Open )
                                                                    year_group = 5,
                                                                    age_max_mort = 110, # Open )
                                                                    age_group = 5), # [,)
                   .id = "Sim_id") 
save(asmr_dir_wd, file = "Measures/asmr_dir_wd.RData")


## Calculate ASFR and ASMR for genealogical subset of direct ancestors without duplicates 

# Keep only unique pids for each simulation and create a list of data frames containing opop of ancestors 
ancestors_dir_wod <- ancestors_10 %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all = TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)

# Estimate age-specific fertility rates for the genealogical subset of direct ancestors without duplicates
asfr_dir_wod <- map_dfr(ancestors_dir_wod, ~ estimate_fertility_rates(opop = .x,
                                                                      final_sim_year = 2022, #[Jan-Dec]
                                                                      year_min = 1750, # Closed [
                                                                      year_max = 2020, # Open )
                                                                      year_group = 5, 
                                                                      age_min_fert = 10, # Closed [
                                                                      age_max_fert = 55, # Open )
                                                                      age_group = 5), # [,)
                             .id = "Sim_id") 
save(asfr_dir_wod, file = "Measures/asfr_dir_wod.RData")

# Estimate age-specific mortality rates for the genealogical subset of direct ancestors without duplicates
asmr_dir_wod <- map_dfr(ancestors_dir_wod, ~ estimate_mortality_rates(opop = .x,
                                                                      final_sim_year = 2022, #[Jan-Dec]
                                                                      year_min = 1750, # Closed
                                                                      year_max = 2020, # Open )
                                                                      year_group = 5,
                                                                      age_max_mort = 110, # Open )
                                                                      age_group = 5), # [,)
                   .id = "Sim_id") 
save(asmr_dir_wod, file = "Measures/asmr_dir_wod.RData")

#----------------------------------------------------------------------------------------------------
## Comparison of whole simulation with subsets of direct ancestors, with and without duplicates ----

# ASFR ----

# Load mean ASFR 5x5 from the 10 simulations, calculated on 2_Compare_Input_Output
load("Measures/asfr_whole.RData")
# Load ASFR for the genealogical subset of direct ancestors with duplicates
load("Measures/asfr_dir_wd.RData")
# Load ASMR for the genealogical subset of direct ancestors with duplicates
load("Measures/asmr_dir_wd.RData")

## Calculate the mean of the different simulations and add relevant columns

# Whole SOCSIM simulation
asfr_whole2 <- asfr_whole %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Whole simulation", 
         Rate = "ASFR") 

# Genealogical subset of direct ancestors with duplicates
asfr_dir_wd2 <- asfr_dir_wd %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct ancestors with duplicates",
         Rate = "ASFR") 

# Genealogical subset of direct ancestors without duplicates
asfr_dir_wod2 <- asfr_dir_wod %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct ancestors without duplicates", 
         Rate = "ASFR") 

## Plot ASFR from whole SOCSIM simulation and the genealogical subset of direct ancestors

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2, asfr_dir_wd2, asfr_dir_wod2) %>% 
  filter(year %in% yrs_plot) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset)))+
  geom_line(aes(colour = year, linetype = Dataset), linewidth = 1.2)+ 
  scale_color_manual(values = c("#79B727", "#B72779", "#2779B7")) +
  scale_linetype_manual(values = c("11", "22", "solid")) +
  theme_graphs()
# labs(title = "Age-Specific Fertility Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and subsets of direct ancestors") 
ggsave(file="Graphs/Socsim_Exp1_ASFR.jpeg", width=17, height=9, dpi=200)


# ASMR ----

# Load mean ASMR rates 5x5 from the 10 simulations, calculated on 2_Compare_Input_Output
load("Measures/asmr_whole.RData")
# Load ASFR for the genealogical subset of direct ancestors without duplicates
load("Measures/asfr_dir_wod.RData")
# Load ASMR for the genealogical subset of direct ancestors without duplicates
load("Measures/asmr_dir_wod.RData")

# Whole SOCSIM simulation
asmr_whole2 <- asmr_whole %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Whole simulation",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)
  
# Genealogical subset of direct ancestors with duplicates
asmr_dir_wd2 <- asmr_dir_wd %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Direct ancestors with duplicates",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Genealogical subset of direct ancestors without duplicates
asmr_dir_wod2 <- asmr_dir_wod %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Direct ancestors without duplicates",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)


## Plotting ASMR from whole SOCSIM simulation and some genealogical subsets

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_dir_wd$age)

bind_rows(asmr_whole2, asmr_dir_wd2, asmr_dir_wod2) %>% 
  rename(Year = year) %>% 
  filter(Year %in% yrs_plot) %>%  
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>% 
  ggplot(aes(x = age, y = mx, group = interaction(Year, Dataset)))+
  facet_wrap(~Sex) +
  geom_line(aes(colour = Year, linetype = Dataset), linewidth = 1.2)+ 
  scale_color_manual(values = c("#79B727", "#B72779", "#2779B7")) +
  scale_linetype_manual(values = c("11", "22", "solid")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs()
  #labs(title = "Age-Specific Mortality Rates in Sweden (1751-2022), 
 # retrieved from a SOCSIM simulation and subset of direct ancestors") 
ggsave(file="Graphs/Socsim_Exp1_ASMR.jpeg", width=17, height=9, dpi=200)
#----------------------------------------------------------------------------------------------------
## Final plot combining ASFR and ASMR ----

# Change years to plot only to two periods
yrs_plot <- c("[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_whole2$age)

## Plotting ASFR and ASMR (for females) from whole SOCSIM simulation and some genealogical subsets
By_Age <- 
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_dir_wd2 %>% rename(Estimate = ASFR),
          asfr_dir_wod2 %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx),
            asmr_dir_wd2 %>% rename(Estimate = mx),
            asmr_dir_wod2 %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels), 
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates")) %>%
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(linewidth = 1.3, show.legend = TRUE)+
  geom_point(aes(shape = Dataset), size = 11)+
  scale_color_manual(values = c( "#B72779", "#2779B7"))+
  scale_shape_manual(values = c(15, 19 ,46)) + 
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(x = "Age") +
  guides(shape = guide_legend(order = 1), col = guide_legend(order = 2)) +
  theme(legend.justification = "left")

# Save the plot
By_Age
ggsave(file="Graphs/Final_Socsim_Exp1_ASFR_ASMR.jpeg", width=17, height=9, dpi=200)

## Save as .svg file for poster
# ggsave(file="Graphs/Socsim_Exp1_ASFR_ASMR.svg", device = "svg", units = "in", width=15, height=8, dpi=200) 

#----------------------------------------------------------------------------------------------------
## Summary measures: TFR and e0 ----
# Here, we use the rates by 1 year age group and 1 calendar year

# Total Fertility Rate ----
# Calculate Total Fertility Rate from asfr 1x1

# Estimate age-specific fertility rates 1x1 for the genealogical subset of direct ancestors with duplicates
# We need to use here the modified function that allows the joining of data frames with duplicates
asfr_dir_wd_1 <- map_dfr(ancestors_dir_wd, ~ estimate_fertility_rates_mod(opop = .x,
                                                                          final_sim_year = 2022, #[Jan-Dec]
                                                                          year_min = 1750, # Closed [
                                                                          year_max = 2023, # Open )
                                                                          year_group = 1, 
                                                                          age_min_fert = 10, # Closed [
                                                                          age_max_fert = 55, # Open )
                                                                          age_group = 1), # [,)
                   .id = "Sim_id") 
save(asfr_dir_wd_1, file = "Measures/asfr_dir_wd_1.RData")

# Estimate age-specific fertility rates 1x1 for the genealogical subset of direct ancestors without duplicates
asfr_dir_wod_1 <- map_dfr(ancestors_dir_wod, ~ estimate_fertility_rates(opop = .x,
                                                                final_sim_year = 2022, #[Jan-Dec]
                                                                year_min = 1750, # Closed [
                                                                year_max = 2023, # Open )
                                                                year_group = 1, 
                                                                age_min_fert = 10, # Closed [
                                                                age_max_fert = 55, # Open )
                                                                age_group = 1), # [,)
                    .id = "Sim_id") 
save(asfr_dir_wod_1, file = "Measures/asfr_dir_wod_1.RData")

# Load ASFR 1x1 and calculate TFR for plotting ----

# Load mean ASFR 1x1 from the 10 simulations, calculated on 2_Compare_Input_Output\
load("Measures/asfr_whole_1.RData")
load("Measures/asfr_dir_wd_1.RData")
load("Measures/asfr_dir_wod_1.RData")

# Age breaks of fertility rates. Extract all the unique numbers from the intervals 
age_breaks_fert_1 <- unique(as.numeric(str_extract_all(asfr_whole_1$age, "\\d+", simplify = T)))

# Estimate age_group size
age_group_fert_1 <- unique(diff(age_breaks_fert_1))

# Whole SOCSIM simulation
TFR_whole <- asfr_whole_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "Whole simulation",
         Rate = "TFR", 
         sex = "female")

# Genealogical subset of direct ancestors with duplicates
TFR_dir_wd <- asfr_dir_wd_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "Direct ancestors with duplicates",
         Rate = "TFR", 
         sex = "female") 

# Genealogical subset of direct ancestors without duplicates
TFR_dir_wod <- asfr_dir_wod_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Dataset = "Direct ancestors without duplicates",
         Rate = "TFR", 
         sex = "female")

## Plott TFR from whole SOCSIM simulation and genealogical subsets of direct ancestors with(out) duplicates
bind_rows(TFR_whole, TFR_dir_wd, TFR_dir_wod) %>% 
  filter(Year >= 1751) %>%
  ggplot(aes(x = Year, y = TFR)) +
  geom_line(aes(colour = Dataset, linetype = Dataset), linewidth = 1.3)+ 
  scale_color_manual(values = c("#7A0005", "#38007A" , "#007A75"))+
  scale_linetype_manual(values = c("11", "22", "solid")) +
  theme_graphs() 
# labs(title = "Total Fertility Rate in Sweden (1751-2022), retrieved from a SOCSIM simulation and genealogical subset of direct ancestors") 
ggsave(file="Graphs/Socsim_Exp1_TFR.jpeg", width=17, height=9, dpi=200)

# Life Expectancy at birth ----
# Estimate life expectancy at birth from asmr 1x1

# Estimate age-specific mortality rates for the genealogical subset of direct ancestor with duplicates
asmr_dir_wd_1 <- map_dfr(ancestors_dir_wd, ~ estimate_mortality_rates(opop = .x,
                                                                      final_sim_year = 2022, #[Jan-Dec]
                                                                      year_min = 1750, # Closed
                                                                      year_max = 2023, # Open )
                                                                      year_group = 1,
                                                                      age_max_mort = 110, # Open )
                                                                      age_group = 1), # [,)
                     .id = "Sim_id") 
save(asmr_dir_wd_1, file = "Measures/asmr_dir_wd_1.RData")

# Estimate age-specific mortality rates for the genealogical subset of direct ancestors without duplicates
asmr_dir_wod_1 <- map_dfr(ancestors_dir_wod, ~ estimate_mortality_rates(opop = .x,
                                                                       final_sim_year = 2022, #[Jan-Dec]
                                                                       year_min = 1750, # Closed
                                                                       year_max = 2023, # Open )
                                                                       year_group = 1,
                                                                       age_max_mort = 110, # Open )
                                                                       age_group = 1), # [,)
                    .id = "Sim_id") 
save(asmr_dir_wod_1, file = "Measures/asmr_dir_wod_1.RData")


# Load mean Age-specific mortality rates 1x1 from the 10 simulations, calculated on 2_Compare_Input_Output
load("Measures/asmr_whole_1.RData")
# Compute life table from mean asmr 1x1 of 10 whole SOCSIM simulations
lt_whole <- lt_socsim(asmr_whole_1)
save(lt_whole, file = "Measures/lt_whole.RData")

# Calculate the mean from asmr 1x1 for the genealogical subset of direct ancestors with duplicates
asmr_dir_wd_1 <- asmr_dir_wd_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
# Compute life table for the genealogical subset of direct ancestors with duplicates
lt_dir_wd <- lt_socsim(asmr_dir_wd_1)
save(lt_dir_wd, file = "Measures/lt_dir_wd.RData")

# Calculate the mean from asmr 1x1 for the genealogical subset of direct ancestors without duplicates
asmr_dir_wod_1 <- asmr_dir_wod_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
# Compute life table for the genealogical subset of direct ancestors without duplicates
lt_dir_wod <- lt_socsim(asmr_dir_wod_1)
save(lt_dir_wod, file = "Measures/lt_dir_wod.RData")

# Load and wrangle life tables for plotting ----
load("Measures/lt_whole.RData")
load("Measures/lt_dir_wd.RData")
load("Measures/lt_dir_wod.RData")

# Year breaks. Extract all the unique numbers from the intervals 
load("Measures/asmr_whole_1.RData")
year_breaks_mort_1 <- unique(as.numeric(str_extract_all(asmr_whole_1$year, "\\d+", simplify = T)))

# Year range to filter HMD data
year_range_mort_1 <- min(year_breaks_mort_1):max(year_breaks_mort_1-1)

# Whole SOCSIM simulation
lt_whole2 <- lt_whole %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Whole simulation",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Genealogical subset of direct ancestors with duplicates
lt_dir_wd2 <- lt_dir_wd %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct ancestors with duplicates",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# Genealogical subset of direct ancestors without duplicates
lt_dir_wod2 <- lt_dir_wod %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct ancestors without duplicates",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

bind_rows(lt_whole2, lt_dir_wd2, lt_dir_wod2) %>% 
  filter(Age == 0 ) %>%
  mutate(Sex = ifelse(sex == "female", "Female", "Male")) %>% 
  ggplot(aes(x = Year, y = ex, group = Dataset))+
  geom_line(aes(colour = Dataset, linetype = Dataset), linewidth = 1.3)+
  scale_color_viridis(option = "D", discrete = T, direction = -1)+
  scale_linetype_manual(values = c("11", "22", "solid")) +
  facet_wrap(~Sex) +
  theme_graphs()+
  labs(title = "Life expectancy at birth in Sweden (e0), 1751-2022, a SOCSIM simulation and genealogical subset of direct ancestors",
       y = "e0") 
ggsave(file="Graphs/Socsim_Exp1_e0.jpeg", width=17, height=9, dpi=200)


#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

## TFR and e0 (for females) from whole SOCSIM simulation and genealogical subsets of direct ancestors

yrs_plot2 <- c(1750, 1800, 1850, 1900, 1950, 2000)

Summary <-
bind_rows(TFR_whole %>% rename(Estimate = TFR), 
          TFR_dir_wd %>% rename(Estimate = TFR),
          TFR_dir_wod %>% rename(Estimate = TFR)) %>%  
  bind_rows(lt_whole2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_dir_wd2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_dir_wod2 %>% rename(Estimate = ex) %>% filter(Age == 0)) %>% 
  filter(sex == "female") %>% 
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life Expectancy at Birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life Expectancy at Birth"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), aes(shape = Dataset), size = 11)+
  geom_line(linewidth = 1.2, show.legend = TRUE) +
  scale_color_manual(values = c("#7A0005", "#38007A" , "#007A75"))+
  scale_shape_manual(values = c(15,19,46)) + 
  scale_x_continuous(breaks = yrs_plot2)+
  theme_graphs() + 
  theme(legend.justification = "left")
# labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2022), retrieved from HFD, HMD and 10 SOCSIM simulation outputs") + 

# Save the plot
Summary
ggsave(file="Graphs/Final_Socsim_Exp1_TFR_e0.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Plot combining age-specific rates and summary measures -----

plot_labs1 <- data.frame(Rate = c("Age-Specific Fertility Rates", "Age-Specific Mortality Rates"),
                         x = c(1,2),
                         y = c(0.19, 0.5),
                         labels = c("a","b"))
plot_labs2 <- data.frame(Rate = as.factor(c("Total Fertility Rate", "Life Expectancy at Birth")),
                         x = c(1755, 1755),
                         y = c(5.2, 82),
                         labels = c("c","d"))

By_Age + 
  geom_text(data = plot_labs1, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") + 
  theme(plot.margin = margin(0,0,1,0, "cm")) +
  Summary + 
  geom_text(data = plot_labs2, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") +
  plot_layout(ncol = 1)
ggsave(file="Graphs/Final_Socsim_Exp1_Combined.jpeg", width=18, height=20, dpi=200)