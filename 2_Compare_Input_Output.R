#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Genealogies - Compare Input and Output Rates
# U:/SOCSIM/SOCSIM_Genealogies/2_Compare_Input_Output.R

## Compare input and output age-specific rates from 10 SOCSIM microsimulations for Sweden (1751-2022)
# as well as summary measures such as TFR, e0

## To run the following code, it is necessary to have already run the simulations and read the .opop file
# c.f. script 1_Run_Simulations.R

# Created by Liliana Calderon on 18-01-2022
# Last modified by Liliana Calderon on 20-07-2023

# NB: Some functions are adapted from external code specified under each section.
#----------------------------------------------------------------------------------------------------
## General settings and functions -----
# Prevent scientific notation (useful for the rate calculation)
options(scipen=999999)

## Load packages 
library(tidyverse)
library(ggh4x)  # To facet scales-
library(HMDHFDplus)
library(patchwork) # To combine ggplots
library(rsocsim) # Functions to estimate rates
library(svglite) # To save svg files
library(viridis)

## Load theme for the graphs and to convert SOCSIM time
source("Functions_Graphs.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with asmr calculated with rsocsim::estimate_mortality_rates()
source("Functions_Life_Table.R")

# Load saved list with opop from simulations, generated in 1_Run_Simulations.R
load("sims_opop.RData")
# Load saved list with omar from simulations, generated in 1_Run_Simulations.R
load("sims_omar.RData")

## Type inside the quotes HFD and HMD credentials (username and password) for the new website
# This information is necessary to run the comparisons below

# HFD credentials
HFD_username <- "Type_here_HFD_username"
HFD_password <- "Type_here_HFD_password"

# HMD credentials
HMD_username <- "Type_here_HMD_username"
HMD_password <- "Type_here_HMD_password"

#------------------------------------------------------------------------------------------------------
# Age-Specific Fertility and Mortality rates, 5x5 ----

# Create a sub-folder called "Measures" to save the output measures if it does not exist.
ifelse(!dir.exists("Measures"), dir.create("Measures"), FALSE)

# Retrieve age-specific fertility rates
asfr_10 <- map_dfr(sims_opop, ~ estimate_fertility_rates(opop = .x,
                                                         final_sim_year = 2022, #[Jan-Dec]
                                                         year_min = 1750, # Closed [
                                                         year_max = 2020, # Open )
                                                         year_group = 5, 
                                                         age_min_fert = 10, # Closed [
                                                         age_max_fert = 55, # Open )
                                                         age_group = 5), # [,)
                .id = "Sim_id") 
save(asfr_10, file = "Measures/asfr_10.RData")

# Retrieve age-specific mortality rates
asmr_10 <- map_dfr(sims_opop, ~ estimate_mortality_rates(opop = .x,
                                                         final_sim_year = 2022, #[Jan-Dec]
                                                         year_min = 1750, # Closed
                                                         year_max = 2020, # Open )
                                                         year_group = 5,
                                                         age_max_mort = 110, # Open )
                                                         age_group = 5), # [,)
                .id = "Sim_id") 
save(asmr_10, file = "Measures/asmr_10.RData")

## Plots ASFR and ASMR

# Load ASFR and ASMR from the 10 simulations
load("Measures/asfr_10.RData")
load("Measures/asmr_10.RData")

# Create a sub-folder called "Graphs" to save the plots if it does not exist.
ifelse(!dir.exists("Graphs"), dir.create("Graphs"), FALSE)

## ASFR and ASMR from SOCSIM for women
# Age groups for fertility and mortality must be the same to plot the figure below

# Choose years to plot (in intervals).
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_10$age)

bind_rows(asfr_10 %>%
            mutate(rate = "ASFR",                   
                   sex = "female"),
          asmr_10 %>% 
            mutate(rate = "ASMR") %>% 
            filter(sex == "female")) %>% 
  mutate(age = factor(as.character(age), levels = age_levels)) %>% 
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(socsim !=0 & !is.infinite(socsim) & !is.nan(socsim)) %>% 
  filter(year %in% yrs_plot) %>% 
  ggplot(aes(x = age, y = socsim, group = interaction(year, Sim_id), colour = year)) +
  geom_line(linewidth = 1) +
  facet_wrap(. ~ rate, scales = "free") + 
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10"))) +
  theme_graphs() +
  scale_color_viridis(option = "D", discrete = T, direction = -1) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "Age", y = "Estimate")
  # labs(title = "Age-Specific Fertility and Mortality Rates for Women retrieved from 10 SOCSIM outputs")
ggsave(file="Graphs/SOCSIM_10_ASFR_ASMR.jpeg", width=17, height=9, dpi=200)


#----------------------------------------------------------------------------------------------------
## Comparison with HFC/HFD and HMD data (used as input) ----

# ASFR ----

# Get HFD ASFR for Sweden (by single year of age and 5 calendar year), 1751-1890
HFC <- read_csv(paste0("https://www.fertilitydata.org/File/GetFile/Country/SWE/SWE_ASFRstand_TOT.txt"),
                col_names = T, show_col_types = F) %>% 
  filter(RefCode %in% "SWE_02" & Year2 <= 1890)

# Get HFD ASFR for Sweden (by single year of age and calendar year), 1891-2022
HFD <- readHFDweb(CNTRY = "SWE",
                  item = "asfrRR",
                  username = HFD_username,
                  password = HFD_password)

# Repeat rates of year groups for each calendar year
HFC <- HFC  %>% 
  select(Year1, Age, ASFR) %>% 
  mutate(ASFR = if_else(ASFR == ".", "0", ASFR), 
         ASFR = as.double(ASFR), 
         Year2 = Year1 + 1, 
         Year3 = Year1 + 2, 
         Year4 = Year1 + 3,
         Year5 = Year1 + 4) %>% 
  select(Year1, Year2, Year3, Year4, Year5, Age, ASFR) %>%
  pivot_longer(cols = c(Year1:Year5), names_to = "Delete", values_to = "Year") %>% 
  select(Year, Age, ASFR) %>% 
  arrange(Year, Age)

# Extract year and age breaks used in the get_asfr_socsim() to apply the same values to HFD data

# Year breaks. Extract all the unique numbers from the intervals. 
year_breaks_fert <- unique(as.numeric(str_extract_all(asfr_10$year, "\\d+", simplify = T)))

# Year range to filter HFD data
year_range_fert <- min(year_breaks_fert):max(year_breaks_fert-1)

# Age breaks of fertility rates. Extract all the unique numbers from the intervals 
age_breaks_fert <- unique(as.numeric(str_extract_all(asfr_10$age, "\\d+", simplify = T)))


# Wrangle HFC and HFD data
HFCD0 <- bind_rows(HFC, HFD) %>% 
  filter(Year %in% year_range_fert) %>% 
  select(-OpenInterval) %>% 
  mutate(year = cut(Year, breaks = year_breaks_fert, 
                       include.lowest = F, right = F, ordered_results = T),
         age = cut(Age, breaks = age_breaks_fert, 
                     include.lowest = F, right = F, ordered_results = T)) %>% 
  filter(!is.na(age)) %>% 
  group_by(year, age) %>%
  summarise(ASFR = mean(ASFR)) %>%
  ungroup() %>%
  mutate(Source = "HFC/HFD", 
         Sim_id = "0", # To plot multiple simulations
         Rate = "ASFR")

# Wrangle SOCSIM data
SocsimF0 <- asfr_10 %>% 
  rename(ASFR = socsim) %>% 
  mutate(Source = "SOCSIM",
         Rate = "ASFR")

## Plot ASFR from HFD vs SOCSIM   

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(HFCD0, SocsimF0) %>%
  filter(year %in% yrs_plot) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Sim_id)))+
  geom_line(aes(colour = year, linetype = Source, alpha = Source), linewidth = 1.2)+
  scale_color_viridis(option = "D", discrete = T, direction = 1) +
  scale_linetype_manual(values = c("HFC/HFD" = "solid", "SOCSIM" = "dotted")) +
  scale_alpha_discrete(guide="none", range = c(1, 0.4))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(title = "Age-Specific Fertility rates in Sweden (1751-2022), retrieved from HFD and 10 SOCSIM simulation outputs") 
ggsave(file="Graphs/HFD_SOCSIM_10_ASFR.jpeg", width=17, height=9, dpi=200)


# ASMR ----
## We compare the mid-year ASMR with central death rates from HMD life tables 1x1
# Rates will be aggregated as defined in the get_asmr_socsim()

# Get female life tables from HMD, 1x1
ltf <- readHMDweb(CNTRY = "SWE",
                  item = "fltper_1x1",
                  username = HMD_username,
                  password = HMD_password)

# Get male life tables from HMD, 1x1
ltm <- readHMDweb(CNTRY = "SWE",
                  item = "mltper_1x1",
                  username = HMD_username,
                  password = HMD_password)

# Extract year and age breaks used in the get_asmr_socsim() to apply the same values to HMD data

# Year breaks. Extract all the unique numbers from the intervals 
year_breaks_mort <- unique(as.numeric(str_extract_all(asmr_10$year, "\\d+", simplify = T)))

# Year range to filter HMD data
year_range_mort <- min(year_breaks_mort):max(year_breaks_mort-1)

# Age breaks of mortality rates. Extract all the unique numbers from the intervals 
age_breaks_mort <- unique(as.numeric(str_extract_all(asmr_10$age, "\\d+", simplify = T)))

# Wrangle HMD life tables
HMD <- ltf %>%
  select(Year, Age, mx) %>% 
  mutate(Sex = "Female") %>% 
  bind_rows(ltm %>% 
              select(Year, Age, mx) %>%  
              mutate(Sex = "Male")) %>% 
  filter(Year %in% year_range_mort) %>% 
  mutate(year = cut(Year, breaks = year_breaks_mort, 
                    include.lowest = F, right = F, ordered_results = T),
         age = cut(Age, breaks = age_breaks_mort, 
                   include.lowest = F, right = F, ordered_results = T)) %>% 
  filter(!is.na(age)) %>% 
  group_by(year, Sex, age) %>% 
    summarise(mx = mean(mx)) %>%
    ungroup() %>%
    mutate(Source = "HMD",
           Sim_id = "0", 
           Rate = "ASMR")

# Wrangle SOCSIM data
SocsimM <- asmr_10 %>% 
  rename(mx = socsim) %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Source = "SOCSIM",
         Rate = "ASMR") %>% 
  select(year, Sex, age,  mx, Source, Sim_id, Rate)


## Plot ASMR from HMD vs SOCSIM   

# Same years to plot than above (in intervals). Change if necessary
# yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(HMD, SocsimM) %>% 
  filter(year %in% yrs_plot) %>% 
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>% 
  ggplot(aes(x = age, y = mx, group = interaction(year, Sim_id))) +
  facet_wrap(~Sex) +
  geom_line(aes(colour = year, linetype = Source, alpha = Source), linewidth = 1.3)+
  scale_y_log10() +
  scale_color_viridis(option = "D", discrete = T, direction = -1) +
  scale_linetype_manual(values = c("HMD" = "solid", "SOCSIM" = "dotted")) +
  scale_alpha_discrete(guide="none", range = c(1, 0.4))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs()+
  labs(title = "Age-Specific Mortality Rates (log-scale) in Sweden (1751-2022), HMD and 10 Socsim simulations (without NaNs and Inf values)")
ggsave(file="Graphs/HMD_SOCSIM_10_log_NA.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Final plot combining ASFR and ASMR ----

# Change years to plot only to two periods
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(SocsimM$age)

## Plotting ASFR and ASMR (for females) from HFD/HMD vs SOCSIM 
By_Age <- 
bind_rows(HFCD0 %>% rename(Estimate = ASFR), 
            SocsimF0 %>% rename(Estimate = ASFR)) %>% 
    mutate(Sex = "Female") %>%   
    bind_rows(HMD %>% rename(Estimate = mx),
              SocsimM %>% rename(Estimate = mx)) %>% 
    # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
    filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
    filter(Sex == "Female") %>% 
    mutate(Year = year,
           age = factor(as.character(age), levels = age_levels), 
           transp = ifelse(Source == "SOCSIM", "0", "1"), 
           Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", "Age-Specific Mortality Rates")) %>%
    filter(Year %in% yrs_plot) %>% 
    ggplot(aes(x = age, y = Estimate, group = interaction(Year, Sim_id)))+
    facet_wrap(. ~ Rate, scales = "free") + 
    geom_line(aes(colour = Year, alpha = transp), linewidth = 1.5)+
    scale_color_manual(values = c("#79B727", "#B72779", "#2779B7"))+
    scale_alpha_discrete(guide="none", range = c(0.2, 1))+
    facetted_pos_scales(y = list("Age-Specific Fertility Rates" = scale_y_continuous(),
                                 "Age-Specific Mortality Rates" =  scale_y_continuous(trans = "log10")))+
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = "Age") +
    # labs(title = "Age-Specific Fertility and Mortality rates in Sweden (1751-2022), retrieved from HFD, HMD and 10 SOCSIM simulation outputs") + 
    theme_graphs()

# Save the plot
By_Age
ggsave(file="Graphs/Final_Socsim_HFD_HMD1.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Summary measures: TFR and e0 ----
# Here, we use the socsim rates by 1 year age group and 1 calendar year

# Total Fertility Rate ----

## Retrieve age-specific fertility rates, by 1 year age group and 1 calendar year
asfr_10_1 <- map_dfr(sims_opop, ~ estimate_fertility_rates(opop = .x,
                                                          final_sim_year = 2022 , #[Jan-Dec]
                                                          year_min = 1750, # Closed [
                                                          year_max = 2023, # Open )
                                                          year_group = 1, 
                                                          age_min_fert = 10, # Closed [
                                                          age_max_fert = 55, # Open )
                                                          age_group = 1), # [,)
                   .id = "Sim_id") 
save(asfr_10_1, file = "Measures/asfr_10_1.RData")
load("Measures/asfr_10_1.RData")

# Year breaks. Extract all the unique numbers from the intervals. 
year_breaks_fert_1 <- unique(as.numeric(str_extract_all(asfr_10_1$year, "\\d+", simplify = T)))

# Year range to filter HFD data
year_range_fert_1 <- min(year_breaks_fert_1):max(year_breaks_fert_1-1)

# Age breaks of fertility rates. Extract all the unique numbers from the intervals 
age_breaks_fert_1 <- unique(as.numeric(str_extract_all(asfr_10_1$age, "\\d+", simplify = T)))

# Retrieve age_group size
age_group_fert_1 <- unique(diff(age_breaks_fert_1))


# Calculate TFR from HFC and HFD
HFCD1 <- bind_rows(HFC, HFD) %>% 
  filter(Year %in% year_range_fert_1) %>% 
  select(-OpenInterval) %>% 
  group_by(Year) %>% 
  summarise(TFR = sum(ASFR, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Source = "HFC/HFD", 
         Sim_id = "0") # To plot multiple simulations

# Calculate TFR from SOCSIM
SocsimF1 <-  asfr_10_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert_1) %>%
  ungroup() %>% 
  mutate(Source = "SOCSIM") 

## Plot TFR from HFD vs SOCSIM 
bind_rows(HFCD1, SocsimF1) %>%
  mutate(transp = ifelse(Source == "SOCSIM", "0", "1")) %>% 
  ggplot(aes(x = Year, y = TFR, group = interaction(Source, Sim_id))) +
  geom_line(aes(colour = Source, alpha = transp), linewidth = 1.3)+
  scale_color_manual(values = c("#007A75", "#CA650D"))+
  scale_alpha_discrete(guide = "none", range = c(0.2, 1))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(title = "Total Fertility rates in Sweden (1751-2022), retrieved from HFD and 10 Socsim simulation outputs") 
ggsave(file="Graphs/HFD_SOCSIM_10_TFR.jpeg", width=17, height=9, dpi=200)


# Life Expectancy at birth ----
# Calculate life expectancy at birth 1x1 for the 10 SOCSIM simulations

# Retrieve age-specific mortality rates, by 1 year age group and 1 calendar year
asmr_10_1 <- map_dfr(sims_opop, ~ estimate_mortality_rates(opop = .x,
                                                          final_sim_year = 2022, #[Jan-Dec]
                                                          year_min = 1750, # Closed
                                                          year_max = 2023, # Open )
                                                          year_group = 1,
                                                          age_max_mort = 110, # Open )
                                                          age_group = 1), # [,)
                     .id = "Sim_id") 
save(asmr_10_1, file = "Measures/asmr_10_1.RData")
load("Measures/asmr_10_1.RData")

# Compute life table from SOCSIM's asmr 1x1 for each simulation
lt_10 <- lt_socsim_sims(asmr_socsim_sims = asmr_10_1)
save(lt_10, file = "Measures/lt_10.RData")
load("Measures/lt_10.RData")

## Compare with ex at age 0 for Sweden in HMD

# Get female life tables from HMD
ltf <- readHMDweb(CNTRY = "SWE",
                  item = "fltper_1x1",
                  username = HMD_username,
                  password = HMD_password)

# Get male life tables from HMD
ltm <- readHMDweb(CNTRY = "SWE",
                  item = "mltper_1x1",
                  username = HMD_username,
                  password = HMD_password)


# Year breaks. Extract all the unique numbers from the intervals 
year_breaks_mort_1 <- unique(as.numeric(str_extract_all(asmr_10_1$year, "\\d+", simplify = T)))

# Year range to filter HMD data
year_range_mort_1 <- min(year_breaks_mort_1):max(year_breaks_mort_1-1)

# Wrangle HMD life tables 
HMD_lt <- ltf %>%
  select(Year, Age, ex) %>% 
  mutate(Sex = "female") %>% 
  bind_rows(ltm %>% 
              select(Year, Age, ex) %>%  
              mutate(Sex = "male")) %>% 
  mutate(Sex = factor(Sex, levels = c("male", "female")), 
         Source = "HMD",
         Sim_id = "0")

# Wrangle SOCSIM life tables
SOCSIM_lt <- lt_10 %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Source = "SOCSIM") %>% 
  select(Year, Age, ex, Sex = sex, Source, Sim_id)

bind_rows(HMD_lt, SOCSIM_lt) %>% 
  filter(Age == 0 & Year %in% year_range_mort_1) %>%
  mutate(transp = ifelse(Source == "SOCSIM", "0", "1"),
         Sex = ifelse(Sex == "female", "Female", "Male")) %>% 
  ggplot(aes(x = Year, y = ex, group = interaction(Source, Sim_id)))+
  geom_line(aes(colour = Source, alpha = transp), linewidth = 1.3)+
  scale_color_manual(values = c("#0C0B7F", "#007A75"))+
  scale_alpha_discrete(guide = "none", range = c(0.2, 1))+
  facet_wrap(~Sex) +
  theme_graphs()+
  labs(title = "Life Expectancy at Birth in Sweden (e0), 1751-2022, retrieved from HMD and 10 SOCSIM simulation outputs",
       y = "e0") 
ggsave(file="Graphs/HMD_SOCSIM_10_e0.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

## Plotting TFR and e0 (for females) from HFD/HMD vs SOCSIM 
Summary <- 
bind_rows(HFCD1 %>% rename(Estimate = TFR) %>%  mutate(Rate = "TFR"),
          SocsimF1 %>% rename(Estimate = TFR) %>% mutate(Rate = "TFR")) %>% 
  mutate(Sex = "female") %>%   
  bind_rows(HMD_lt %>% filter(Age == 0) %>% rename(Estimate = ex) %>% mutate(Rate = "e0"),
            SOCSIM_lt %>% filter(Age == 0) %>% rename(Estimate = ex) %>% mutate(Rate = "e0")) %>% 
  filter(Sex == "female") %>%
  mutate(transp = ifelse(Source == "SOCSIM", "0", "1"),
         Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life Expectancy at Birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life Expectancy at Birth"))) %>% 
  ggplot(aes(x = Year, y = Estimate, group = interaction(Source, Sim_id)))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(aes(colour = Source, alpha = transp), linewidth = 1.2) +
  scale_color_manual(values = c("#CA650D", "#0C0B7F", "#007A75")) +
  scale_alpha_discrete(guide="none", range = c(0.1, 1)) +
  scale_x_continuous(breaks = c(1750, 1800, 1850, 1900, 1950, 2000))+
  theme_graphs()
  # labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2022), retrieved from HFD, HMD and 10 SOCSIM simulation outputs") + 
# Save the plot
Summary
ggsave(file="Graphs/Final_Socsim_HFD_HMD2.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
#### Plot combining age-specific rates and summary measures -----

plot_labs1 <- data.frame(Rate = c("Age-Specific Fertility Rates", "Age-Specific Mortality Rates"),
                        x = c(1,2),
                        y = c(0.22, 0.5),
                        labels = c("a","b"))
plot_labs2 <- data.frame(Rate = as.factor(c("Total Fertility Rate", "Life Expectancy at Birth")),
                        x = c(1755, 1755),
                        y = c(5.3, 81),
                        labels = c("c","d"))

By_Age + 
  geom_text(data = plot_labs1, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") + 
  theme(plot.margin = margin(0,0,1,0, "cm")) +
  Summary + 
  geom_text(data = plot_labs2, mapping = aes(x = x, y = y, label = labels), inherit.aes = F, 
            size = 15, family="serif") +
  plot_layout(ncol = 1)

ggsave(file="Graphs/Final_Socsim_HFD_HMD_Combined.jpeg", width=18, height=20, dpi=200)

#----------------------------------------------------------------------------------------------------
## Calculate the mean measures of the different simulations ----

# Age-specific fertility rates 5x5
load("Measures/asfr_10.RData")
asfr_whole <- asfr_10 %>%
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
save(asfr_whole, file = "Measures/asfr_whole.RData")

# Age-specific mortality rates 5x5
load("Measures/asmr_10.RData")
asmr_whole <- asmr_10 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
save(asmr_whole, file = "Measures/asmr_whole.RData")

# Age-specific fertility rates 1x1
load("Measures/asfr_10_1.RData")
asfr_whole_1 <- asfr_10_1 %>%
  group_by(year, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
save(asfr_whole_1, file = "Measures/asfr_whole_1.RData")

# Age-specific mortality rates 1x1
load("Measures/asmr_10_1.RData")
asmr_whole_1 <- asmr_10_1 %>%
  group_by(year, sex, age) %>% 
  summarise(socsim = mean(socsim, na.rm = T)) %>% 
  ungroup()
save(asmr_whole_1, file = "Measures/asmr_whole_1.RData")