#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Genealogies - Compare Input and Output Rates
# U:/SOCSIM/SOCSIM_Genealogies/2_Compare_Input_Output.R

## Compare input and output age-specific rates from 10 SOCSIM microsimulations for Sweden (1751-2021)
# The script also includes some comparisons of summary measures like TFR, e0, SRB, IMR

## To run the following code, it is necessary to have already run the simulations and read the .opop file
# c.f. script 1_Run_Simulations.R

# Created by Liliana Calderon on 18-01-2022
# Last modified by Liliana Calderon on 14-03-2023

# NB: Some functions are adapted from external code specified under each section.

#----------------------------------------------------------------------------------------------------
## Recovering the input mortality and fertility rates -----

# Prevent scientific notation (useful for the rate calculation)
options(scipen=999999)

# Load packages 
library(tidyverse)

# Load saved list with 10 simulations, generated in 1_Run_Simulations.R
load("sims_opop.RData")

## Load functions to recover input age-specific fertility and mortality rates 
# These are a modified version of scripts `02_extract_rates.R` and `functions.R` written by Diego Alburez-Gutierrez,
# available on https://github.com/alburezg/socsim_example/tree/main/R 
source("Functions_Retrieve_Rates.R")

# Retrieve age-specific fertility rates
asfr_10 <- map_dfr(sims_opop, ~ get_asfr_socsim(df = .x,
                                             final_sim_year = 2021 , #[Jan-Dec]
                                             year_min = 1750, # Closed [
                                             year_max = 2020, # Open )
                                             year_group = 5, 
                                             age_min_fert = 10, # Closed [
                                             age_max_fert = 55, # Open )
                                             age_group = 5), # [,)
                .id = "Sim_id") 
save(asfr_10, file = "asfr_10.RData")

# Retrieve age-specific mortality rates
asmr_10 <- map_dfr(sims_opop, ~ get_asmr_socsim(df = .x,
                                             final_sim_year = 2021, #[Jan-Dec]
                                             year_min = 1750, # Closed
                                             year_max = 2020, # Open )
                                             year_group = 5,
                                             age_max_mort = 110, # Open )
                                            age_group = 5), # [,)
                .id = "Sim_id") 
save(asmr_10, file = "asmr_10.RData")

#----------------------------------------------------------------------------------------------------
## Plot the results ----
library(viridis)
library(ggh4x) # For facetted_pos_scales
library(svglite) # To save svg files

## Theme for the graphs
source("Functions_Graphs.R")

# Load ASFR and ASMR from the 10 simulations
load("asfr_10.RData")
load("asmr_10.RData")

# Create a sub-folder called "Graphs" to save the plots if it does not exist.
ifelse(!dir.exists("Graphs"), dir.create("Graphs"), FALSE)

## ASFR and ASMR from SOCSIM for women ----
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
ggsave(file="Graphs/SOCSIM_10_ASFR_ASMR.jpeg", width=17, height=9, dpi=400)


#----------------------------------------------------------------------------------------------------
## Comparison with HFD and HMD data (used as input) ----

## Load HMDHFD package to access the data
# install.packages("HMDHFDplus")
library(HMDHFDplus)

#### Age-Specific Fertility rates ----

# HFD credentials
HFD_username <- "Type_here_HFD_username"
HFD_password <- "Type_here_HFD_password"

# HFD ASFR for Sweden (by calendar year and single year of age)
HFD <- readHFDweb(CNTRY = "SWE",
                  item = "asfrRR",
                  username = HFD_username,
                  password = HFD_password)

# Extract year and age breaks used in the get_asfr_socsim() to apply the same values to HFD data

# Year breaks. Extract all the unique numbers from the intervals. 
year_breaks_fert <- unique(as.numeric(str_extract_all(asfr_10$year, "\\d+", simplify = T)))

# Year range to filter HFD data
year_range_fert <- min(year_breaks_fert):max(year_breaks_fert-1)

# Age breaks of fertility rates. Extract all the unique numbers from the intervals 
age_breaks_fert <- unique(as.numeric(str_extract_all(asfr_10$age, "\\d+", simplify = T)))

# Wrangle HFD data
HFD0 <- HFD %>% 
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
  mutate(Source = "HFD", 
         Sim_id = "0", # To plot multiple simulations
         Rate = "ASFR")

# Wrangle SOCSIM data
SocsimF0 <- asfr_10 %>% 
  rename(ASFR = socsim) %>% 
  mutate(Source = "SOCSIM",
         Rate = "ASFR")

## Plot ASFR from HFD vs SOCSIM   

# Same years to plot than above (in intervals). Change if necessary
# yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(HFD0, SocsimF0) %>%
  # HFD rates for [1890,1895) are used for[1800,1805) in the simulation
  mutate(Year = case_when(Source == "HFD" & year == "[1890,1895)" ~ "[1800,1805)", 
                          TRUE ~ year)) %>%
  filter(Year %in% yrs_plot) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(Year, Sim_id)))+
  geom_line(aes(colour = Year, linetype = Source, alpha = Source), linewidth = 1.2)+
  scale_color_viridis(option = "D", discrete = T, direction = 1) +
  scale_linetype_manual(values = c("HFD" = "solid", "SOCSIM" = "dotted")) +
  scale_alpha_discrete(guide="none", range = c(1, 0.4))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(title = "Age-Specific Fertility rates in Sweden (1751-2021), retrieved from HFD and 10 SOCSIM simulation outputs") 
ggsave(file="Graphs/HFD_SOCSIM_10_ASFR.jpeg", width=17, height=9, dpi=400)


####  Age-Specific Mortality rates ----
## We compare the mid-year ASMR with central death rates from HMD life tables 1x1
# Rates will be aggregated as defined in the get_asmr_socsim()

# HMD credentials
HMD_username <- "Type_here_HMD_username"
HMD_password <- "Type_here_HMD_password"

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
  labs(title = "Age-Specific Mortality Rates (log-scale) in Sweden (1751-2021), HMD and 10 Socsim simulations (without NaNs and Inf values)")
ggsave(file="Graphs/HMD_SOCSIM_10_log_NA.jpeg", width=17, height=9, dpi=400)


## Final plot combining ASFR and ASMR ----

# Same years to plot than above (in intervals). Change if necessary
# yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(SocsimM$age)

## Plotting ASFR and ASMR (for females) from HFD/HMD vs SOCSIM 
bind_rows(HFD0 %>% rename(Estimate = ASFR), 
            SocsimF0 %>% rename(Estimate = ASFR)) %>% 
    mutate(Sex = "Female") %>%   
    bind_rows(HMD %>% rename(Estimate = mx),
              SocsimM %>% rename(Estimate = mx)) %>% 
    # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
    filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
    filter(Sex == "Female") %>% 
    # Hack to set HFD rates for 1890 as 1800
    mutate(Year = case_when(Source == "HFD"& Rate == "ASFR" & year == "[1890,1895)" ~ "[1800,1805)", 
                            TRUE ~ year),
           age = factor(as.character(age), levels = age_levels), 
           transp = ifelse(Source == "SOCSIM", "0", "1")) %>%
    filter(Year %in% yrs_plot) %>% 
    ggplot(aes(x = age, y = Estimate, group = interaction(Year, Sim_id)))+
    facet_wrap(. ~ Rate, scales = "free") + 
    geom_line(aes(colour = Year, linetype = Source, alpha = transp), linewidth = 1.5)+
    scale_color_manual(values = c("#35B779", "#31688E", "#440154"))+
    scale_linetype_manual(values = c("HFD" = "solid", "HMD" = "solid","SOCSIM" = "dotted")) +
    scale_alpha_discrete(guide="none", range = c(0.5, 1))+
    facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                                 ASMR =  scale_y_continuous(trans = "log10")))+
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = "Age") +
    # labs(title = "Age-Specific Fertility and Mortality rates in Sweden (1751-2021), retrieved from HFD, HMD and 10 SOCSIM simulation outputs") + 
    theme_graphs() 

# Save the plot
ggsave(file="Graphs/Final_Socsim_HFD_HMD1.jpeg", width=17, height=9, dpi=400)

#----------------------------------------------------------------------------------------------------
#### Summary measures: TFR and e0 ----
# Here, we use the rates by single calendar year and year of age

#### Total Fertility Rate

## Retrieve age-specific fertility rates, by single calendar year and 1 year age group
asfr_10_1 <- map_dfr(sims_opop, ~ get_asfr_socsim(df = .x,
                                                final_sim_year = 2021 , #[Jan-Dec]
                                                year_min = 1750, # Closed [
                                                year_max = 2020, # Open )
                                                year_group = 1, 
                                                age_min_fert = 10, # Closed [
                                                age_max_fert = 55, # Open )
                                                age_group = 1), # [,)
                   .id = "Sim_id") 
save(asfr_10_1, file = "asfr_10_1.RData")

load("asfr_10_1.RData")

# Year breaks. Extract all the unique numbers from the intervals. 
year_breaks_fert <- unique(as.numeric(str_extract_all(asfr_10_1$year, "\\d+", simplify = T)))

# Year range to filter HFD data
year_range_fert <- min(year_breaks_fert):max(year_breaks_fert-1)

# Age breaks of fertility rates. Extract all the unique numbers from the intervals 
age_breaks_fert <- unique(as.numeric(str_extract_all(asfr_10_1$age, "\\d+", simplify = T)))

# Retrieve age_group size
age_group_fert <- unique(diff(age_breaks_fert))

# Calculate TFR from HFD
HFD1 <- HFD %>% 
  filter(Year %in% year_range_fert) %>% 
  select(-OpenInterval) %>% 
  group_by(Year) %>% 
  summarise(TFR = sum(ASFR)) %>% 
  ungroup() %>% 
  mutate(Source = "HFD", 
         Sim_id = "0") 

# Calculate TFR from SOCSIM
SocsimF1 <-  asfr_10_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year, Sim_id) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Source = "SOCSIM") 

## Plot TFR from HFD vs SOCSIM 
bind_rows(HFD1, SocsimF1) %>% 
  mutate(transp = ifelse(Source == "SOCSIM", "0", "1")) %>% 
  ggplot(aes(x = Year, y = TFR, group = interaction(Source, Sim_id))) +
  geom_line(aes(colour = Source, alpha = transp), linewidth = 1.3)+
  scale_color_manual(values = c("#30734A", "#CA650D"))+
  scale_alpha_discrete(guide = "none", range = c(0.2, 1))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(title = "Total Fertility rates in Sweden (1751-2021), retrieved from HFD and 10 Socsim simulation outputs") 
ggsave(file="Graphs/HFD_SOCSIM_10_TFR.jpeg", width=17, height=9, dpi=400)


## Life Expectancy at birth ----
# Calculate life expectancy at birth 1x1 for the 10 SOCSIM simulations

# Retrieve age-specific mortality rates, by single calendar year and 1 year age group
asmr_10_1 <- map_dfr(sims_opop, ~ get_asmr_socsim(df = .x,
                                                  final_sim_year = 2021, #[Jan-Dec]
                                                  year_min = 1750, # Closed
                                                  year_max = 2020, # Open )
                                                  year_group = 1,
                                                  age_max_mort = 110, # Open )
                                                  age_group = 1), # [,)
                     .id = "Sim_id") 
save(asmr_10_1, file = "asmr_10_1.RData")

load("asmr_10_1.RData")

# This code was inspired by Tim Riffe's BSSD2021Module2 code to calculate life tables
# https://github.com/timriffe/BSSD2021Module2/blob/master/02_tuesday/02_tuesday.Rmd 

lt_sim <- asmr_10_1 %>% 
  mutate(Age = as.numeric(str_extract(age, "\\d+"))) %>% 
  rename(mx = socsim) %>% 
  group_by(Sim_id, year, sex) %>% 
  mutate(n = ifelse(Age == max(Age), 1, lead(Age)-Age), # n = 1 for Open Age Interval as in rate files
         rn = row_number()) %>% 
  # Filter data frame until maximum row with mx greater than 0
  filter(between(rn, 1, max(which(mx > 0)))) %>% 
  mutate(mx = case_when(is.nan(mx) ~ 0,
                        is.infinite(mx) ~ 1,  # Set Infinite mx(x/0) to 1
                        TRUE ~ mx),
         # Use a0 formulas from HMD (Modified version from Andreev and Kingkade, 2015)
         ax = case_when(Age == 0 & sex == "female" & mx < 0.01724 ~ 0.14903 - 2.05527 * mx,
                        Age == 0 & sex == "female" & mx < 0.06891 ~ 0.04667 + 3.88089 * mx,
                        Age == 0 & sex == "female" & mx >= 0.06891 ~ 0.31411,
                        Age == 0 & sex == "male" & mx < 0.02300 ~ 0.14929 - 1.99545 * mx,
                        Age == 0 & sex == "male" & mx < 0.08307 ~ 0.02832 + 3.26021 * mx,
                        Age == 0 & sex == "male" & mx >= 0.08307 ~ 0.29915,
                        Age == max(Age) ~ 1/mx, # With this transformation, there can be very high ax
                        TRUE ~ n/2),
         qx = (mx * n) / (1 + (n - ax) * mx), 
         qx = case_when(Age == max(Age) ~ 1, 
                        qx > 1 ~ 1, TRUE ~ qx)) %>%  
  # Filter data frame until minimum row with qx == 1, to avoid having more than one qx = 1
  filter(between(rn, 1, min(which(qx ==1)))) %>%
  mutate(px = 1 - qx,
         lx = 1e5 * c(1, cumprod(px[-n()])),
         dx = lx * qx,
         Lx = n * lx - (n - ax) * dx, 
         Tx = Lx %>% rev() %>% cumsum() %>% rev(),
         ex = Tx / lx)  %>% 
  ungroup() 
save(lt_sim, file = "lt_sim.RData")


## Compare with ex at age 0 for Sweden in HDM
library(HMDHFDplus)

# HMD credentials
HMD_username <- "Type_here_HMD_username"
HMD_password <- "Type_here_HMD_password"

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
SOCSIM_lt <- lt_sim %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Source = "SOCSIM") %>% 
  select(Year, Age, ex, Sex = sex, Source, Sim_id)

bind_rows(HMD_lt, SOCSIM_lt) %>% 
  filter(Age == 0 & Year %in% year_range_mort_1) %>%
  mutate(transp = ifelse(Source == "SOCSIM", "0", "1"),
         Sex = ifelse(Sex == "female", "Female", "Male")) %>% 
  ggplot(aes(x = Year, y = ex, group = interaction(Source, Sim_id)))+
  geom_line(aes(colour = Source, alpha = transp), linewidth = 1.3)+
  scale_color_manual(values = c("#0C0B7F", "#30734A"))+
  scale_alpha_discrete(guide = "none", range = c(0.2, 1))+
  facet_wrap(~Sex) +
  theme_graphs()+
  labs(title = "Life expectancy at birth in Sweden (e0), 1751-2020, retrieved from HMD and 10 SOCSIM simulation outputs",
       y = "e0") 
ggsave(file="Graphs/HMD_SOCSIM_10_e0.jpeg", width=17, height=9, dpi=400)

#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

## Ploting TFR and e0 (for females) from HFD/HMD vs SOCSIM 
bind_rows(HFD1 %>% 
            rename(Estimate = TFR) %>% 
            mutate(Rate = "TFR"),
          SocsimF1 %>% 
            rename(Estimate = TFR) %>% 
            mutate(Rate = "TFR")) %>% 
  mutate(Sex = "female") %>%   
  bind_rows(HMD_lt %>% 
              filter(Age == 0) %>% 
              rename(Estimate = ex) %>% 
              mutate(Rate = "e0"),
            SOCSIM_lt %>% 
              filter(Age == 0) %>% 
              rename(Estimate = ex) %>% 
              mutate(Rate = "e0")) %>% 
  filter(Sex == "female") %>% 
  mutate(transp = ifelse(Source == "SOCSIM", "0", "1"),
         Rate = factor(Rate, levels = c("TFR", "e0"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = interaction(Source, Sim_id)))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(aes(colour = Source, alpha = transp), linewidth = 1.2) +
  scale_color_manual(values = c("#CA650D", "#0C0B7F", "#30734A")) +
  scale_alpha_discrete(guide="none", range = c(0.1, 1)) +
  theme_graphs()
  # labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2020), retrieved from HFD, HMD and 10 SOCSIM simulation outputs") + 
   
# Save the plot
ggsave(file="Graphs/Final_Socsim_HFD_HMD2.jpeg", width=17, height=9, dpi=400)

#----------------------------------------------------------------------------------------------------
## Births and Deaths counts by calendar year -----

last_month <- map_dbl(sims_opop, ~ .x %>% 
                       pull(dob) %>% 
                       max()) %>%
  unique()

## If not set in the Global Environment
final_sim_year <- 2021 #[Jan-Dec]
year_min <- 1750 # Closed [
year_max <- 2020 # Open )

# Year range
year_range <- year_min:(year_max-1)

# Birth counts by year
Births <- map_dfr(sims_opop, ~ .x %>% 
                    mutate(Year = asYr(dob, last_month, final_sim_year)) %>% 
                    filter(Year %in% year_range) %>% 
                    count(Year) %>%
                    mutate(Year = factor(Year, levels = year_range)) %>% 
                    complete(Year, fill = list(n = 0))  %>%
                    mutate(Year = as.numeric(as.character(Year))) %>% 
                    rename(Births = n),
                  .id = "Sim_id")


# Death counts
Deaths <- map_dfr(sims_opop, ~ .x %>% 
                      filter(dod != 0) %>% 
                      mutate(Year = asYr(dod, last_month, final_sim_year)) %>% 
                      filter(Year %in% year_range) %>% 
                      count(Year) %>%
                      mutate(Year = factor(Year, levels = year_range)) %>% 
                      complete(Year, fill = list(n = 0))  %>%
                      mutate(Year = as.numeric(as.character(Year))) %>% 
                      rename(Deaths = n),
                    .id = "Sim_id")

# Plotting birth and death counts together
full_join(Births, Deaths,  by = c("Sim_id", "Year")) %>% 
  filter(Year >= 1751 & Sim_id %in% c("1", "2", "3")) %>%
  mutate(Sim_id = factor(Sim_id, levels = 1:10)) %>% 
  pivot_longer(Births:Deaths, names_to = "Event", values_to = "Count") %>% 
  ggplot(aes(x = Year, y = Count, group = interaction(Sim_id), color = Sim_id))+
  facet_wrap(. ~ Event) + 
  geom_line(linewidth = 1)+
  scale_color_viridis(option = "C", discrete = T, direction = -1) +
  theme_graphs()
ggsave(file="Graphs/Socsim_Births_Deaths.jpeg", width=17, height=9, dpi=400)

#----------------------------------------------------------------------------------------------------
## Sex Ratio at Birth and Infant Mortality Rate
# We use here the asYr() function from the Functions_Retrieve_Rates.R
## Check if this should be included in a function

## If not set in the Global Environment
final_sim_year <- 2021 #[Jan-Dec]
year_min <- 1750 # Closed [
year_max <- 2020 # Open )

# Year range
year_range <- year_min:(year_max-1)

## Sex Ratio at Birth 

# Birth counts by year and sex
Births_Sex <- map_dfr(sims_opop, ~ .x %>% 
                        mutate(Year = asYr(dob, last_month, final_sim_year), 
                               Sex = ifelse(fem == 1, "Female", "Male")) %>% 
                        filter(Year %in% year_range) %>% 
                        count(Year, Sex) %>%
                        mutate(Year = factor(Year, levels = year_range),
                               Sex = factor(Sex, levels = c("Female", "Male"))) %>% 
                        complete(Year, Sex, fill = list(n = 0))  %>%
                        mutate(Year = as.numeric(as.character(Year)),
                               Sex = as.character(Sex)) %>% 
                        rename(Births = n),
                      .id = "Sim_id")

# Calculate and plot Sex Ratio at Birth (SRB)
SRB <- Births_Sex %>% 
  pivot_wider(names_from = Sex, values_from = Births) %>% 
  mutate(SRB = Male/Female, 
         Sim_id = factor(Sim_id, levels = 1:10)) 

SRB %>% 
  filter(Year >= 1751) %>% 
  ggplot(aes(x = Year, y = SRB, group = Sim_id, color = Sim_id))+
  geom_line(linewidth = 1)+
  scale_color_viridis(option = "C", discrete = T, direction = -1) +
  theme_graphs()


#### Infant Mortality Rate, both sexes

# Birth counts by year
Births <- map_dfr(sims_opop, ~ .x %>% 
                    mutate(Year = asYr(dob, last_month, final_sim_year)) %>% 
                    filter(Year %in% year_range) %>% 
                    count(Year) %>%
                    mutate(Year = factor(Year, levels = year_range)) %>% 
                    complete(Year, fill = list(n = 0))  %>%
                    mutate(Year = as.numeric(as.character(Year))) %>% 
                    rename(Births = n),
                  .id = "Sim_id")


# Death counts below age 1 (0-11 months)
Deaths_0 <- map_dfr(sims_opop, ~ .x %>% 
                      filter(dod != 0) %>% 
                      mutate(age_death_months = dod-dob,
                             Year = asYr(dod, last_month, final_sim_year)) %>% 
                      filter(age_death_months < 12 & Year %in% year_range) %>% 
                      count(Year) %>%
                      mutate(Year = factor(Year, levels = year_range)) %>% 
                      complete(Year, fill = list(n = 0))  %>%
                      mutate(Year = as.numeric(as.character(Year))) %>% 
                     rename(Deaths = n),
                    .id = "Sim_id")

# Calculate and Plot Infant Mortality Rate (IMR)
IMR <- full_join(Births, Deaths_0, 
                 by = c("Sim_id", "Year")) %>% 
  mutate(IMR = Deaths/Births, 
         Sim_id = factor(Sim_id, levels = 1:10)) 

IMR %>% 
  filter(Year >= 1751) %>% 
  ggplot(aes(x = Year, y = IMR, group = Sim_id, color = Sim_id))+
  geom_line(linewidth = 1)+
  scale_color_viridis(option = "C", discrete = T, direction = -1) +
  theme_graphs()


# Plot SRB and IMR together

full_join(SRB %>% select(Sim_id, Year, SRB),
          IMR %>% select(Sim_id, Year, IMR),
          by = c("Sim_id", "Year")) %>% 
  filter(Year >= 1751 & Sim_id %in% c("1", "2", "3")) %>%
  mutate(Sim_id = factor(Sim_id, levels = 1:10)) %>% 
  pivot_longer(SRB:IMR, names_to = "Measure", values_to = "Value") %>% 
  mutate(Measure = factor(Measure, levels = c("SRB", "IMR"))) %>%
    ggplot(aes(x = Year, y = Value, group = Sim_id, color = Sim_id))+
  facet_wrap(. ~ Measure, scales = "free") + 
  geom_line(linewidth = 1)+
  scale_color_viridis(option = "C", discrete = T, direction = -1) +
  facetted_pos_scales(y = list(SRB = scale_y_continuous(limits=c(0.75, 1.25)),
                               IMR =  scale_y_continuous())) +
  theme_graphs()
ggsave(file="Graphs/Socsim_SRB_IMR.jpeg", width=17, height=9, dpi=400)

#----------------------------------------------------------------------------------------------------
#### Infant Mortality Rate by sex. 
## To compare with ASMR at age 0

# Death counts by sex below age 1 (0-11 months)
Deaths_Sex_0 <- map_dfr(sims_opop, ~ .x %>% 
                      filter(dod != 0) %>% 
                      mutate(age_death_months = dod-dob,
                             Year = asYr(dod, last_month, final_sim_year),
                             Sex = ifelse(fem == 1, "Female", "Male")) %>% 
                      filter(age_death_months < 12 & Year %in% year_range) %>% 
                      count(Year, Sex) %>%
                      mutate(Year = factor(Year, levels = year_range),
                             Sex = factor(Sex, levels = c("Female", "Male"))) %>% 
                      complete(Year, Sex, fill = list(n = 0))  %>%
                      mutate(Year = as.numeric(as.character(Year)),
                             Sex = as.character(Sex)) %>% 
                      rename(Deaths = n),
                    .id = "Sim_id")

# Calculate and Plot Infant Mortality Rate (IMR by sex)
IMR_sex <- full_join(Births_Sex, Deaths_Sex_0, 
          by = c("Sim_id", "Year", "Sex")) %>% 
  mutate(IMR = Deaths/Births, 
         Sim_id = factor(Sim_id, levels = 1:10)) 

IMR_sex %>% 
  filter(Year >= 1751) %>% 
  ggplot(aes(x = Year, y = IMR, group = Sim_id, color = Sim_id))+
  facet_wrap(. ~ Sex) + 
  geom_line(linewidth = 1)+
  scale_color_viridis(option = "C", discrete = T, direction = -1) +
  theme_graphs()

# Compare with death rates at age 0 from asmr
load("asmr_10_1.RData")

# Rates from asmr (using mid-year pop aged 0) are a bit higher than IMR (using births as denominator)
asmr_10_1 %>% 
  filter(age == "[0,1)") %>% 
  mutate(Sim_id = factor(Sim_id, levels = 1:10),
         Year = as.numeric(str_extract(year, "\\d+")),
         Sex = ifelse(sex =="female", "Female", "Male"),
         Source = "Soc1") %>% 
  select(Sim_id, Year, Sex, IMR = socsim, Source) %>% 
  filter(Year >= 1751) %>%
  ggplot(aes(x = Year, y = IMR, group = Sim_id, color = Sim_id))+
  facet_wrap(. ~ Sex) + 
  geom_line(linewidth = 1)+
  scale_color_viridis(option = "C", discrete = T, direction = -1) +
  theme_graphs()