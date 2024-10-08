#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Genealogies - Trace and compare subset of direct ancestors of a given ego(s) 
# U:/SOCSIM/SOCSIM_Genealogies/3_Compare_Ancestors.R

## Trace direct ancestors of a given ego(s) from a SOCSIM microsimulation for Sweden (1751-2072)
# and compare demographic measures from the whole simulation and the genealogical subsets

# Created on 23-09-2022
# Last modified on 03-09-2024

## NB: To run this code, it is necessary to have already run the script 1_Run_Simulations.R
#------------------------------------------------------------------------------------------------------
## General settings and functions ----

# Prevent scientific notation (useful for the rate calculation)
options(scipen=999999)

## Load packages 
library(tidyverse)
library(ggh4x)  # To facet scales-
library(patchwork) # To combine ggplots

## Load function to get direct ancestors
source("Functions/Functions_Ancestors.R")

## Load functions to estimate age-specific fertility and mortality rates 
# These are a slightly modified version of the functions in the rsocsim package 
# that allow to handle the intentional duplicates in the data (in fertility rates)
# and retrieve the last month from max(dod) instead of dob. 
# Important here as max(dob) is not the last simulated month
# So, the function in rsocsim would assign incorrectly the last month of the simulation
# and hence calculate wrongly the years of birth and death
# Also, the direct ancestors offspring (e.g., here, only ego, but siblings/aunts/uncles later) 
# are only counted in the numerator but not in the denominator of fertility rates
source("Functions/Functions_Fertility_Rates_Mod.R")
source("Functions/Functions_Mortality_Rates_Mod.R")

## Load theme for the graphs and to convert SOCSIM time
source("Functions/Functions_Graphs.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with asmr calculated with estimate_mortality_rates()
source("Functions/Functions_Life_Table.R")
#------------------------------------------------------------------------------------------------------
## Trace direct ancestors of people alive in 2023 ----

# Load saved list with opop from 10 simulations, generated in 1_Run_Simulations.R
load("sims_opop.RData")
# Load saved list with omar from 10 simulations, generated in 1_Run_Simulations.R
load("sims_omar.RData")
# Create a sub-folder called "Subsets" to save the opop subsets used for each experiment
ifelse(!dir.exists("Subsets"), dir.create("Subsets"), FALSE)

# Function to get a 10% sample of people alive at the end of 2022

sample_egos <- function(opop = opop, final_sim_year, percentage) {
  egos <- opop %>% 
    mutate(last_month = max(dob),
           birth_year = asYr(dob, last_month, final_sim_year),
           death_year = asYr(dod, last_month, final_sim_year)) %>%  
    # People born before 2022 and alive by the end of 2022
    filter(birth_year <= final_sim_year-50 & death_year > final_sim_year-50) %>% 
    pull(pid)
  sample_size <- round(length(egos)*percentage/100)
  egos_samp <-  sample(egos, sample_size, replace = F)
  return(egos_samp)
}

# Get 10% sample for each of the 10 simulations
egos_samp_10 <- map(sims_opop, ~ sample_egos(opop = .x,
                                             final_sim_year = 2072, 
                                             percentage = 10)) 
save(egos_samp_10, file = "Subsets/egos_samp_10.RData")


start <- Sys.time()
## Retrieve the ancestors of each simulation sample of egos alive in 2023
ancestors_10 <- map2_dfr(egos_samp_10, sims_opop,
                           ~ retrieve_ancestors(egos = .x, opop = .y), 
                           .id = "Sim_id") 
end <- Sys.time()
print(end-start) # Time difference of 2.588096 mins
# Save the data frame with the ancestors of 10 simulations samples
save(ancestors_10, file = "Subsets/ancestors_10.RData")

#----------------------------------------------------------------------------------------------------
## Check range of births of egos and female ancestors

last_month <- ancestors_10 %>% 
  pull(dod) %>% 
  (max)

final_sim_year <- 2072

## Check year range of ego's dob
ancestors_10 %>% 
  filter(kin_type == "ego") %>% 
  pull(dob) %>% 
  range() %>% 
  asYr(., last_month, final_sim_year) # 1917-2022

## Check year range of female ancestors's dob
ancestors_10 %>% 
  filter(kin_type != "ego" & fem ==1) %>% 
  pull(dob) %>% 
  range() %>% 
  asYr(., last_month, final_sim_year) # 1602-1999
# So, there might be some female ancestors for births in the last years, but aged at least 22

# The max(dob) is not the last simulated month
# So, the function in rsocsim will assign incorrectly the last month of the simulation
# and hence calculate wrongly the years of birth and death
# We must use the modified rate functions (using dod) to ensure that the last month correctly retrieved
# With max(dod), in all subsets and simulations, the last month will be properly identified

#----------------------------------------------------------------------------------------------------
## Age-Specific Fertility and Mortality rates, 5x5  -----
# Estimate and compare the rates derived from the whole simulation with those from the genealogical subset
# of direct ancestors both with and without duplicates

# Load the data frame with the ancestors of the 10 simulations samples of egos alive by the end of 2022
load("Subsets/ancestors_10.RData")

#  Estimate ASFR and ASMR for genealogical subset of direct ancestors of population alive in 01-01-2023 with duplicates 

# Create a list of data frames containing opop of ancestors for each simulation
ancestors_dir_wd <- ancestors_10 %>%
  split(.$Sim_id)

## Check maximum dod to check if the last month is correctly retrieved (with the modified version)
ancestors_dir_wd %>% 
  map_dbl(.f = ~ max(.x$dod)) %>% 
  unique() %>% 
  range() %>% 
  asYr(., last_month, final_sim_year)

# Estimate age-specific fertility rates for the genealogical subset of Direct Ancestors (with duplicates)
# We need to use here the modified function that allows the joining of data frames with duplicates
# And extract last_month from max(dod) instead of max(dob)
asfr_dir_wd <- map_dfr(ancestors_dir_wd, ~ estimate_fertility_rates_mod(opop = .x,
                                                                        final_sim_year = 2072, #[Jan-Dec]
                                                                        year_min = 1750, # Closed [
                                                                        year_max = 2020, # Open )
                                                                        year_group = 5, 
                                                                        age_min_fert = 10, # Closed [
                                                                        age_max_fert = 55, # Open )
                                                                        age_group = 5), # [,)
                       .id = "Sim_id") 
save(asfr_dir_wd, file = "Measures/asfr_dir_wd.RData")

# Estimate age-specific mortality rates for the genealogical subset of direct ancestor with duplicates
# We need to use here the modified function that extracts last_month from max(dod) instead of max(dob)
asmr_dir_wd <- map_dfr(ancestors_dir_wd, ~ estimate_mortality_rates_mod(opop = .x,
                                                                        final_sim_year = 2072, #[Jan-Dec]
                                                                        year_min = 1750, # Closed
                                                                        year_max = 2020, # Open )
                                                                        year_group = 5,
                                                                        age_max_mort = 110, # Open )
                                                                        age_group = 5), # [,)
                   .id = "Sim_id") 
save(asmr_dir_wd, file = "Measures/asmr_dir_wd.RData")


## Calculate ASFR and ASMR for genealogical subset of Direct Ancestors (without duplicates) 

# Keep only unique pids for each simulation and create a list of data frames containing opop of ancestors 
ancestors_dir_wod <- ancestors_10 %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all = TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)

# Estimate age-specific fertility rates for the genealogical subset of Direct Ancestors (without duplicates)
asfr_dir_wod <- map_dfr(ancestors_dir_wod, ~ estimate_fertility_rates_mod(opop = .x,
                                                                          final_sim_year = 2072, #[Jan-Dec]
                                                                          year_min = 1750, # Closed [
                                                                          year_max = 2020, # Open )
                                                                          year_group = 5, 
                                                                          age_min_fert = 10, # Closed [
                                                                          age_max_fert = 55, # Open )
                                                                          age_group = 5), # [,)
                             .id = "Sim_id") 
save(asfr_dir_wod, file = "Measures/asfr_dir_wod.RData")

# Estimate age-specific mortality rates for the genealogical subset of Direct Ancestors (without duplicates)
asmr_dir_wod <- map_dfr(ancestors_dir_wod, ~ estimate_mortality_rates_mod(opop = .x,
                                                                          final_sim_year = 2072, #[Jan-Dec]
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

# Load ASFR 5x5 from the 10 simulations
load("Measures/asfr_10.RData")
# Load ASFR for the genealogical subset of Direct Ancestors (with duplicates)
load("Measures/asfr_dir_wd.RData")
# Load ASFR for the genealogical subset of Direct Ancestors (without duplicates)
load("Measures/asfr_dir_wod.RData")

## Calculate the mean of the different simulations and add relevant columns

# Whole SOCSIM simulations
asfr_whole2 <- asfr_10 %>%
  group_by(year, age) %>% 
  summarise(ASFR = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Dataset = "Whole Simulation", 
         Rate = "ASFR") 

# Genealogical subset of Direct Ancestors (with duplicates)
asfr_dir_wd2 <- asfr_dir_wd %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct Ancestors (with duplicates)",
         Rate = "ASFR") 

# Genealogical subset of Direct Ancestors (without duplicates)
asfr_dir_wod2 <- asfr_dir_wod %>% 
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Direct Ancestors (without duplicates)", 
         Rate = "ASFR") 

## Plot ASFR from whole SOCSIM simulation and the genealogical subset of direct ancestors

# Same years to plot than above (in intervals). 
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2, asfr_dir_wd2, asfr_dir_wod2) %>% 
  filter(year %in% yrs_plot) %>% 
  mutate(ASFR = ifelse(is.na(ASFR), 0, ASFR)) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset)))+
  geom_line(aes(colour = year, linetype = Dataset), linewidth = 1.2)+ 
  scale_color_manual(values = c("#79B727", "#2779B7", "#B72779")) +
  scale_linetype_manual(values = c("11", "22", "solid")) +
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp1_ASFR_0.jpeg", width=17, height=9, dpi=300)


# ASMR ----

# Load ASMR 5x5 from the 10 simulations
load("Measures/asmr_10.RData")
# Load ASMR for the genealogical subset of Direct Ancestors (without duplicates)
load("Measures/asmr_dir_wod.RData")
# Load ASMR for the genealogical subset of Direct Ancestors (with duplicates)
load("Measures/asmr_dir_wd.RData")

# Whole SOCSIM simulations
asmr_whole2 <- asmr_10 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Whole Simulation",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Genealogical subset of Direct Ancestors (with duplicates)
asmr_dir_wd2 <- asmr_dir_wd %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Direct Ancestors (with duplicates)",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# Genealogical subset of Direct Ancestors (without duplicates)
asmr_dir_wod2 <- asmr_dir_wod %>% 
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Direct Ancestors (without duplicates)",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)


## Plotting ASMR from whole SOCSIM simulation and some genealogical subsets

# Same years to plot than above (in intervals). 
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
  geom_line(aes(colour = Year, linetype = Dataset), linewidth = 1.3)+ 
  scale_linetype_manual(values = c("11", "22", "solid")) +
  scale_color_manual(values = c("#79B727", "#2779B7", "#B72779")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs() +
  labs(x = "Age")
ggsave(file="Graphs/Socsim_Exp1_ASMR.jpeg", width=17, height=9, dpi=300)

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

## Plotting ASFR and ASMR (for females) from whole SOCSIM simulation and some genealogical subsets
By_Age_Exp1 <- 
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_dir_wd2 %>% rename(Estimate = ASFR),
          asfr_dir_wod2 %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx),
            asmr_dir_wd2 %>% rename(Estimate = mx),
            asmr_dir_wod2 %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot1) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels), 
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates")) %>%
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(linewidth = 1.3, show.legend = T)+
  geom_point(data = . %>% filter(age %in% age_plot), 
             aes(shape = Dataset), size = 11) +
  scale_color_manual(values = c("#2779B7"))+ 
  scale_shape_manual(values = c(15, 19 ,46)) + 
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(x = "Age") +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  guides(colour = "none")
By_Age_Exp1

## Plot ASFR and ASMR (for females), with three years for appendix

# Choose three years to plot
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels for fertility to fix them across all plots
age_levels_asfr <- levels(asfr_whole2$age)

bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_dir_wd2 %>% rename(Estimate = ASFR),
          asfr_dir_wod2 %>% rename(Estimate = ASFR)) %>% 
  complete(year, age, Dataset, Rate, fill = list(Estimate = NA)) %>% 
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
  facet_wrap(Year ~ Rate, nrow = 3, ncol = 2, scales = "free") + 
  geom_line(linewidth = 1.3, show.legend = TRUE)+
  geom_point(data = . %>% filter(age %in% age_plot), 
           aes(shape = Dataset), size = 11) +
  scale_color_manual(values = c("#79B727", "#2779B7", "#B72779"))+ 
  scale_shape_manual(values = c(15, 19 ,46)) + 
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10"),
                               ASFR = scale_y_continuous(breaks = y_breaks_asfr),
                               ASMR =  scale_y_continuous(breaks = y_breaks_asmr, trans = "log10"),
                               ASFR = scale_y_continuous(breaks = y_breaks_asfr, limits = c(0, 0.2)),
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
        legend.text = element_text(size = 18)) #or 17
ggsave(file="Final_Graphs/App_Socsim_Exp1_ASFR_ASMR.jpeg", width=18, height=25, dpi=300)

#----------------------------------------------------------------------------------------------------
# Figure for EPC presentation

# Choose three years to plot
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 
age_plot <- c("[0,1)", "[1,5)", "[10,15)", "[20,25)", "[30,35)", "[40,45)", "[50,55)",  "[60,65)", 
              "[70,75)", "[80,85)", "[90,95)", "[100,105)") 

bind_rows(asmr_whole2, asmr_dir_wd2, asmr_dir_wod2) %>% 
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
  scale_shape_manual(values = c(15, 19 ,46)) + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs() +
  labs(x = "Age") +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))+
  guides(colour = "none")

ggsave(file="Graphs/Socsim_Exp1_ASMR_years.jpeg", width=24, height=9, dpi=300)

#----------------------------------------------------------------------------------------------------
## Summary measures: TFR and e0 ----
# Here, we use the rates by 1 year age group and 1 calendar year

# Total Fertility Rate ----
# Estimate Total Fertility Rate from asfr 1x1 ----

# Estimate age-specific fertility rates 1x1 for the genealogical subset of Direct Ancestors (with duplicates)
# We need to use here the modified function that allows the joining of data frames with duplicates
asfr_dir_wd_1 <- map_dfr(ancestors_dir_wd, ~ estimate_fertility_rates_mod(opop = .x,
                                                                          final_sim_year = 2072, #[Jan-Dec]
                                                                          year_min = 1750, # Closed [
                                                                          year_max = 2023, # Open )
                                                                          year_group = 1, 
                                                                          age_min_fert = 10, # Closed [
                                                                          age_max_fert = 55, # Open )
                                                                          age_group = 1), # [,)
                   .id = "Sim_id") 
save(asfr_dir_wd_1, file = "Measures/asfr_dir_wd_1.RData")

# Estimate age-specific fertility rates 1x1 for the genealogical subset of Direct Ancestors (without duplicates)
asfr_dir_wod_1 <- map_dfr(ancestors_dir_wod, ~ estimate_fertility_rates_mod(opop = .x,
                                                                            final_sim_year = 2072, #[Jan-Dec]
                                                                            year_min = 1750, # Closed [
                                                                            year_max = 2023, # Open )
                                                                            year_group = 1, 
                                                                            age_min_fert = 10, # Closed [
                                                                            age_max_fert = 55, # Open )
                                                                            age_group = 1), # [,)
                    .id = "Sim_id") 
save(asfr_dir_wod_1, file = "Measures/asfr_dir_wod_1.RData")

# Load ASFR 1x1 and calculate TFR for plotting ----

# Load asfr 1x1 from the 10 simulations
load("Measures/asfr_10_1.RData")
# Load asfr 1x1 from the subset of direct ancestors with duplicates
load("Measures/asfr_dir_wd_1.RData")
# Load asfr 1x1 from the subset of direct ancestors without duplicates
load("Measures/asfr_dir_wod_1.RData")

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

# Direct Ancestors (with duplicates)
TFR_dir_wd <- asfr_dir_wd_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  # Some ages can have infinite (N_Births/0_Pop) and NaN (0_Births/0_Pop) values
  filter(!is.infinite(socsim)) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim, na.rm = T)) %>%
  ungroup() %>% 
  mutate(Dataset = "Direct Ancestors (with duplicates)",
         Rate = "TFR", 
         sex = "female") 

# Direct Ancestors (without duplicates)
TFR_dir_wod <- asfr_dir_wod_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  # Some ages can have infinite (N_Births/0_Pop) and NaN (0_Births/0_Pop) values
  filter(!is.infinite(socsim)) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim, na.rm = T)) %>%
  ungroup() %>% 
  mutate(Dataset = "Direct Ancestors (without duplicates)",
         Rate = "TFR", 
         sex = "female") 

## Plot TFR from whole SOCSIM simulation and genealogical subsets of direct ancestors with(out) duplicates

yrs_plot2 <- c(1750, 1800, 1850, 1900, 1950, 2000)

bind_rows(TFR_whole, TFR_dir_wd, TFR_dir_wod) %>% 
  group_by(Year, Dataset) %>% 
  summarise(TFR = mean(TFR, na.rm = T)) %>% 
  ungroup() %>% 
  filter(Year >= 1751 & TFR != 0) %>% 
  ggplot(aes(x = Year, y = TFR, colour = Dataset)) +
  geom_point(data = . %>% filter(Year %in% yrs_plot2), aes(shape = Dataset), size = 11)+
  geom_line(linewidth = 1.3, show.legend = TRUE) +
  scale_color_manual(values = c("#00057A", "#7A7500" , "#007A75"))+
  scale_shape_manual(values = c(15,19,46)) + 
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp1_TFR.jpeg", width=17, height=9, dpi=300)

# Summary measure of error in TFR ----

# Difference in means
DiM_TFR_Exp1 <- bind_rows(TFR_whole, TFR_dir_wd,  TFR_dir_wod) %>%
  filter(Year > 1750) %>% 
  # There can be rates of 0 after 1995
  #filter(TFR != 0) %>% 
  group_by(Year, Dataset) %>% 
  summarise(TFR = mean(TFR, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = Year, names_from = "Dataset", values_from = "TFR") %>% 
  pivot_longer(cols = 2:3, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation` , 
         Relative_Error = (Error/`Whole Simulation`)*100,
         Type = "DiM") %>% 
  select(-c(Genealogy,`Whole Simulation`)) 

# Mean of differences
MoD_TFR_Exp1 <- bind_rows(TFR_whole, TFR_dir_wd,  TFR_dir_wod) %>%
  filter(Year > 1750) %>% 
  # There can be rates of 0 after 1995
  #filter(TFR != 0) %>% 
  pivot_wider(id_cols = c(Year, Sim_id), names_from = "Dataset", values_from = "TFR") %>% 
  pivot_longer(cols = 4:5, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation`, 
         Relative_Error = (Error/`Whole Simulation`)*100) %>% 
  group_by(Year, Dataset) %>% 
  summarise(Error = mean(Error, na.rm = T), 
          Relative_Error = mean(Relative_Error, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Type = "MoD")

# Bind both error measures and save the data frame
error_TFR_exp1 <- bind_rows(DiM_TFR_Exp1, MoD_TFR_Exp1) 
save(error_TFR_exp1, file = "Measures/error_TFR_exp1.RData")

# Absolute Error  
error_TFR_exp1 %>% 
  ggplot(aes(x = Year, y = Error, group = Dataset, colour = Type)) +
  facet_wrap(. ~ Dataset)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp1_TFR_Error.jpeg", width=17, height=9, dpi=300)

# Relative Error
error_TFR_exp1 %>% 
ggplot(aes(x = Year, y = Relative_Error, group = Dataset, colour = Type)) +
  facet_wrap(. ~ Dataset)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp1_TFR_Rel_Error.jpeg", width=17, height=9, dpi=300)


# Check mean values of bias in TFR before fertility decline
error_TFR_exp1 %>% 
  filter(!is.na(Error)) %>% # NA after 1997
  filter(Year < 1901) %>% 
  filter(Dataset == "Direct Ancestors (with duplicates)") %>% # -3.453683 -76.91802
  # filter(Dataset == "Direct Ancestors (without duplicates)") %>% # -1.899224 -42.37596%
 # pull(Error) %>% 
  pull(Relative_Error) %>%
  mean()


# Check mean values of bias in TFR 
error_TFR_exp1 %>% 
  filter(!is.na(Error)) %>% # NA after 1997
  #filter(Dataset == "Direct Ancestors (with duplicates)") %>% # -2.447281 -65.78808
  filter(Dataset == "Direct Ancestors (without duplicates)") %>% # -1.467246 -41.71113
  pull(Error) %>% 
 # pull(Relative_Error) %>%
  mean()


# Life Expectancy at birth ----
# Estimate life expectancy at birth from asmr 1x1 for the different genealogical subsets ----

# Estimate age-specific mortality rates for the genealogical subset of direct ancestor with duplicates
asmr_dir_wd_1 <- map_dfr(ancestors_dir_wd, ~ estimate_mortality_rates_mod(opop = .x,
                                                                          final_sim_year = 2072, #[Jan-Dec]
                                                                          year_min = 1750, # Closed
                                                                          year_max = 2023, # Open )
                                                                          year_group = 1,
                                                                          age_max_mort = 110, # Open )
                                                                          age_group = 1), # [,)
                     .id = "Sim_id") 
save(asmr_dir_wd_1, file = "Measures/asmr_dir_wd_1.RData")

# Compute life tables for the genealogical subset of Direct Ancestors (with duplicates)
lt_dir_wd <- lt_socsim_sims(asmr_socsim_sims = asmr_dir_wd_1)
save(lt_dir_wd, file = "Measures/lt_dir_wd.RData")

# Estimate age-specific mortality rates for the genealogical subset of Direct Ancestors (without duplicates)
asmr_dir_wod_1 <- map_dfr(ancestors_dir_wod, ~ estimate_mortality_rates_mod(opop = .x,
                                                                           final_sim_year = 2072, #[Jan-Dec]
                                                                           year_min = 1750, # Closed
                                                                           year_max = 2023, # Open )
                                                                           year_group = 1,
                                                                           age_max_mort = 110, # Open )
                                                                           age_group = 1), # [,)
                    .id = "Sim_id") 
save(asmr_dir_wod_1, file = "Measures/asmr_dir_wod_1.RData")

# Compute life tables for the genealogical subset of Direct Ancestors (without duplicates)
lt_dir_wod <- lt_socsim_sims(asmr_socsim_sims = asmr_dir_wod_1)
save(lt_dir_wod, file = "Measures/lt_dir_wod.RData")

# Load and wrangle life tables for plotting ----

# Load life tables from each whole SOCSIM simulation
load("Measures/lt_10.RData")
# Load life tables from the subset of direct ancestors (with duplicates)
load("Measures/lt_dir_wd.RData")
# Load life tables from the subset of direct ancestors (without duplicates)
load("Measures/lt_dir_wod.RData")

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

# Genealogical subset of Direct Ancestors (with duplicates)
lt_dir_wd2 <- lt_dir_wd %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct Ancestors (with duplicates)",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)

# Genealogical subset of Direct Ancestors (without duplicates)
lt_dir_wod2 <- lt_dir_wod %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Direct Ancestors (without duplicates)",
         Rate = "e0") %>% 
  select(Year, Sim_id, ex, Dataset, Rate, sex, Age)


# Plot the estimates of life expectancy at birth
bind_rows(lt_whole2, lt_dir_wd2, lt_dir_wod2) %>% 
  filter(Age == 0 ) %>%
  mutate(Sex = ifelse(sex == "female", "Female", "Male")) %>% 
  group_by(Year, Sex, Dataset) %>% 
  summarise(ex = mean(ex, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Year, y = ex, group = Dataset))+
  geom_line(aes(colour = Dataset), linewidth = 1.3)+
  scale_color_manual(values = c("#00057A", "#7A7500" , "#007A75"))+
  facet_wrap(~Sex) +
  theme_graphs()+
  labs(y = "Life expectancy at birth")
ggsave(file="Graphs/Socsim_Exp1_e0.jpeg", width=17, height=9, dpi=300)

# Summary measure of error in e0 ----

# Difference in means
DiM_e0_Exp1 <- bind_rows(lt_whole2, lt_dir_wd2, lt_dir_wod2) %>%
  filter(Year > 1750 & Age == 0) %>% 
  group_by(Year, sex, Dataset) %>% 
  summarise(ex = mean(ex, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(Year:sex), names_from = "Dataset", values_from = "ex") %>% 
  pivot_longer(cols = 3:4, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation` , 
         Relative_Error = (Error/`Whole Simulation`)*100,
         Type = "DiM") %>% 
  select(-c(Genealogy,`Whole Simulation`)) 

# Mean of differences
MoD_e0_Exp1 <- bind_rows(lt_whole2, lt_dir_wd2, lt_dir_wod2) %>%
  filter(Year > 1750 & Age == 0) %>% 
  pivot_wider(id_cols = c(Year, sex, Sim_id), names_from = "Dataset", values_from = "ex") %>% 
  pivot_longer(cols = 5:6, names_to = "Dataset", values_to = "Genealogy") %>% 
  mutate(Error = Genealogy - `Whole Simulation`, 
         Relative_Error = (Error/`Whole Simulation`)*100) %>% 
  group_by(Year, sex, Dataset) %>% 
  summarise(Error = mean(Error, na.rm = T), 
          Relative_Error = mean(Relative_Error, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Type = "MoD")

# Bind both error measures and save the data frame
error_e0_exp1 <- bind_rows(DiM_e0_Exp1, MoD_e0_Exp1) 
save(error_e0_exp1, file = "Measures/error_e0_exp1.RData")

# Absolute error
error_e0_exp1 %>% 
  ggplot(aes(x = Year, y = Error, group = Dataset, colour = Type)) +
  facet_grid(Dataset ~ sex)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp1_e0_Error.jpeg", width=17, height=9, dpi=300)

# Relative error
error_e0_exp1 %>% 
  ggplot(aes(x = Year, y = Relative_Error, group = Dataset, colour = Type)) +
  facet_grid(Dataset ~ sex)+
  geom_line(linewidth = 1.3)+
  theme_graphs()
ggsave(file="Graphs/Socsim_Exp1_e0_Rel_Error.jpeg", width=17, height=9, dpi=300)

# Check minimum and maximum values of bias in e0
error_e0_exp1 %>% 
  filter(sex == "female" ) %>% 
  filter(!is.na(Error)) %>% # There can be NAs in error due to missing values
  # filter(Dataset == "Direct Ancestors (with duplicates)") %>% # 0.6294378 34.0463697
  filter(Dataset == "Direct Ancestors (without duplicates)") %>% # 0.7335692 32.9182953
  pull(Error) %>%
  range()

# Check minimum and maximum values of relative bias in e0
error_e0_exp1 %>% 
  filter(sex == "female") %>% 
  filter(!is.na(Error)) %>% # There can be NAs in error due to missing values
  #filter(Dataset == "Direct Ancestors (with duplicates)") %>% # 33.40046
  filter(Dataset == "Direct Ancestors (without duplicates)") %>% # 34.21708
  pull(Relative_Error) %>% 
  mean()

#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

## TFR and e0 (for females) from whole SOCSIM simulation and genealogical subsets of direct ancestors

yrs_plot2 <- c(1750, 1800, 1850, 1900, 1950, 2000)

# Define the same y breaks for all plots
y_breaks_TFR <- c(0:5)
y_breaks_e0 <- c(20, 40, 60, 80)

Summary_Exp1 <-
bind_rows(TFR_whole %>% rename(Estimate = TFR), 
          TFR_dir_wd %>% rename(Estimate = TFR),
          TFR_dir_wod %>% rename(Estimate = TFR)) %>%  
  bind_rows(lt_whole2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_dir_wd2 %>% rename(Estimate = ex) %>% filter(Age == 0),
            lt_dir_wod2 %>% rename(Estimate = ex) %>% filter(Age == 0)) %>% 
  filter(sex == "female") %>% 
  group_by(Year, Dataset, Rate, sex) %>% 
  summarise(Estimate = mean(Estimate, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life Expectancy at Birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life Expectancy at Birth"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), aes(shape = Dataset), size = 11)+
  geom_line(linewidth = 1.2, show.legend = TRUE) +
  scale_color_manual(values = c("#00057A", "#7A7500" , "#007A75"))+
  scale_shape_manual(values = c(15,19,46)) + 
  scale_x_continuous(breaks = yrs_plot2)+
  facetted_pos_scales(y = list("Total Fertility Rate" = scale_y_continuous(breaks = y_breaks_TFR, 
                                                                           limits = c(0, NA)),
                               "Life Expectancy at Birth" =  scale_y_continuous(breaks = y_breaks_e0)))+
  theme_graphs() + 
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

# Save the plot
Summary_Exp1
ggsave(file="Graphs/Socsim_Exp1_TFR_e0.jpeg", width=17, height=9, dpi=300)

#----------------------------------------------------------------------------------------------------
## Plot combining age-specific rates and Summary measures -----

plot_labs1 <- data.frame(Rate = c("Age-Specific Fertility Rates", "Age-Specific Mortality Rates"),
                         x = c(1,2),
                         y = c(0.23, 0.55),
                         labels = c("a","b"))
plot_labs2 <- data.frame(Rate = as.factor(c("Total Fertility Rate", "Life Expectancy at Birth")),
                         x = c(1755, 1755),
                         y = c(5.4, 84),
                         labels = c("c","d"))

By_Age_Exp1 + 
  geom_text(data = plot_labs1, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") + 
  theme(plot.margin = margin(0,0,1,0, "cm")) +
  Summary_Exp1 + 
  geom_text(data = plot_labs2, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") +
  plot_layout(ncol = 1)
ggsave(file="Final_Graphs/Final_Socsim_Exp1_Combined.jpeg", width=18, height=21, dpi=300)