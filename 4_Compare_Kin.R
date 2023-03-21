#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Sweden - Trace and compare a subset of kin of a given ego(s) 
# U:/SOCSIM/SOCSIM_Genealogies/4_Compare_Kin.R

## Trace relevant kin (up to the 4th generation) 
# of a given ego(s) from a SOCSIM microsimulation for Sweden (1751-2021)
# and compare demographic measures from the whole simulation and the genealogical subsets

# Created by Liliana Calderon on 13-04-2022
# Last modified by Liliana Calderon on 17-03-2023

## NB: To run this code, it is necessary to have already run the script 1_Run_Simulations.R

#------------------------------------------------------------------------------------------------------
## General settings ----

# Prevent scientific notation (useful for the rate calculation)
options(scipen=999999)

## Load packages 
library(tidyverse)
library(ggh4x)  # To facet scales-
library(questionr)

## Load functions to recover age-specific fertility and mortality rates 
# and covert SOCSIM time to calendar time
source("Functions_Retrieve_Rates.R")

## Load functions to get kin up to the 4th degree of consanguinity
source("Functions_Kin.R")

## Load theme for the graphs
source("Functions_Graphs.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with the asmr df calculated with the get_asmr_socsim()
source("Functions_Life_Table.R")

# Load read_opop() function to read the .opop file (written by Diego Alburez-Gutierrez)
# Once this is integrated into the rsocsim package the function might be just called
source("read_opop.R")

# Load read_omar() function to read the .omar file (written by Diego Alburez-Gutierrez)
# Once this is integrated into the rsocsim package the function might be just called
source("read_omar.R")

#------------------------------------------------------------------------------------------------------
## Read the output .opop and .omar file ----
## We use only one of the 10 simulations. 

# Name of the supervisory file used for the simulation (if not set in the GlobalEnv)
supfile <- "socsim_SWE.sup"

## Choose simulation seed defined in the previous script
seed <- "13486"

# Path of the simulation results .opop file
path_opop <- paste0("sim_results_",supfile,"_",seed,"_/result.opop")

# Read SOCSIM .opop output, using read_opop function. 
opop <- read_opop(path_opop)

# Path of the simulation results .omar file
path_opop <- paste0("sim_results_",supfile,"_",seed,"_/result.omar")

# Read SOCSIM .omar output, using read_omar function. 
omar <- read_omar(path_opop)

#------------------------------------------------------------------------------------------------------
## Trace kin (up to 4th degree of consanguinity) of people alive in 2022 as a proxy of current genealogists 

# Pids of people alive at the end of the simulation, i.e. dod == 0, 
# who are older than 18 years old on 01-01-2022, i,e. dob 1913-2003 
# egos2022 <- opop %>% 
#   mutate(last_month = max(dob),
#          final_sim_year = 2021, ## Change if necessary
#          Generation = asYr(dob, last_month, final_sim_year)) %>% 
#   filter(dod == 0 & Generation <= final_sim_year-18) %>% 
#   pull(pid)

# Get a sample of 10% of people alive in 2022. 
# sample_size <- round(length(egos2022)/10)
# egos2022_samp <-  sample(egos2022, sample_size, replace = F)
# save(egos2022_samp, file = "egos2022_samp_10.RData")

# Load same sample used to trace the ancestors
load("egos2022_samp_10.RData")

## Map the function to get the ancestors of a sample of individuals alive in 2022 (older than 18 years)
start <- Sys.time()
kin_egos2022_10 <- map_dfr(egos2022_samp, get_kin) %>%
  left_join(select(opop, c(pid, fem, dob, dod, mom, marid, mstat)), by = "pid")
end <- Sys.time()
print(end-start)

# Save the data frame
save(ancestors_egos2022_10, file = "ancestors_egos2022_10.RData")


#----------------------------------------------------------------------------------------------------
# CHECK this later
# Distribution of egos by birth cohort (generation)
bind_rows(kin_alive) %>% 
  filter(pid == ego_id) %>% 
  mutate(cohort = asYr(dob)) %>% 
  pull(cohort) %>% 
  freq()

# Unique pids
bind_rows(kin_alive) %>% 
  pull(pid) %>% 
  unique() %>% 
  length()

# Duplicated pids
bind_rows(kin_alive) %>% 
  pull(pid) %>% 
  duplicated() %>% 
  sum()

# Kin with duplicates for people alive in 2021
kin_alive_wd <- bind_rows(Kin_Alive)

# Kin with duplicates for people alive in 2021
kin_alive_wod <- bind_rows(Kin_Alive) %>% 
  distinct(pid, .keep_all = TRUE)

kin_alive_wod %>% 
  pull(ego_id) %>% 
  unique() %>% 
  length()

#----------------------------------------------------------------------------------------------------
## Recovering the input mortality and fertility rates -----

## To run the following code lines, it is necessary to have already run the script 1_Run_SWE_Example.R
# as it reads Socsim output, attachs the datasets to the search path 
# and defines some necessary constants and functions. 

## Get a glimpse on the data 
kin_alive_wd %>% head()
kin_alive_wod %>% head()

## Adding mom column and calculate years of birth and death using the asYr() function. 

kin_alive_wd <- kin_alive_wd %>%
  left_join(opop %>% 
              select(pid, mom, marid, mstat), 
            by = "pid") %>% 
  mutate(birth_year = asYr(dob),
         death_year = ifelse(dod == 0, NA, asYr(dod))) 

kin_alive_wod <- kin_alive_wod %>%
  left_join(opop %>% 
              select(pid, mom, marid, mstat), 
            by = "pid") %>% 
  mutate(birth_year = asYr(dob),
         death_year = ifelse(dod == 0, NA, asYr(dod))) 

# Define some parameters
y_min <- 1920
y_max <- 2020
y_range <- y_min:y_max
y_breaks <- seq(y_min, y_max, 5)

# Age categories of mortality rates (lower age bounds)
age_breaks_mort <- c(0, 1, seq(5, 100, by = 5))
age_labels_mort <- age_breaks_mort[-length(age_breaks_mort)]

# Age categories of fertility rates (lower age bounds)
age_group_size <-  5
# Reproductive period (Modified to compare with HFD age range -12 to 55+)
min_age <-  10 
max_age <-  55
age_breaks_fert <- seq(min_age, max_age, by = age_group_size)
age_labels_fert <- age_breaks_fert[-length(age_breaks_fert)]


# Load functions to recover input age-specific fertility and mortality rates 
source("Functions_Retrieve_Rates.R")

# Extract age-specific fertility and mortality rates  for the subset with duplicates
asfr_wd <- get_asfr_socsim(df = kin_alive_wd, 
                        y_range = y_range,
                        age_breaks = age_breaks_fert,
                        age_labels = age_labels_fert,
                        sex_keep = "female")

asmr_wd <- get_asmr_socsim(df = kin_alive_wd, 
                        age_breaks = age_breaks_mort, 
                        age_labels = age_labels_mort, 
                        y_breaks = y_breaks,
                        y_range = y_range, 
                        only_women = F)


# Extract age-specific fertility and mortality rates for the subset without duplicates
asfr_wod <- get_asfr_socsim(df = kin_alive_wod, 
                           y_range = y_range,
                           age_breaks = age_breaks_fert,
                           age_labels = age_labels_fert,
                           sex_keep = "female")

asmr_wod <- get_asmr_socsim(df = kin_alive_wod, 
                           age_breaks = age_breaks_mort, 
                           age_labels = age_labels_mort, 
                           y_breaks = y_breaks,
                           y_range = y_range, 
                           only_women = F)

#----------------------------------------------------------------------------------------------------
## Plotting results for subset without duplicates ----
library(viridis)

## ASFR and ASMR without duplicates

bind_rows(asfr_wod %>% 
            mutate(rate = "ASFR", 
                   age = as.numeric(as.character(age))),
          asmr_wod %>% 
            mutate(rate = "ASMR") %>% 
            filter(sex == "female")) %>% 
  mutate(Year = factor(year, 
                       levels = c("[1920,1925]", "(1925,1930]", 
                                  "(1930,1935]", "(1935,1940]", 
                                  "(1940,1945]", "(1945,1950]",
                                  "(1950,1955]", "(1955,1960]", 
                                  "(1960,1965]", "(1965,1970]",
                                  "(1970,1975]", "(1975,1980]", 
                                  "(1980,1985]", "(1985,1990]", 
                                  "(1990,1995]", "(1995,2000]",
                                  "(2000,2005]", "(2005,2010]", 
                                  "(2010,2015]", "(2015,2020]"))) %>% 
  ## Filtering rates = 0 and NAs
  filter(socsim !=0 & !is.na(socsim)) %>% 
  ggplot(aes(x = age, y = socsim, colour = Year)) +
  geom_line(size=1) +
  facet_wrap(. ~ rate, scales = "free") + 
  theme_bw() +
  scale_color_viridis(option = "F", discrete = T) + 
  scale_y_log10()


## Comparison with HFD and HMD data ----

## Load HMDHFD package to access the data
# install.packages("HMDHFDplus")
library(HMDHFDplus)

## Fertility rates ----

# HFD credentials
# HFD_username <- "Type_here_HFD_username"
# HFD_password <- "Type_here_HFD_password"


# HFD ASFR for Sweden (by calendar year and single year of age)
HFD <- readHFDweb(CNTRY = "SWE",
                  item = "asfrRR",
                  username = HFD_username,
                  password = HFD_password)


### Age-specific Fertility Rates

# Wrangling HFD data
HFD0 <- HFD %>% 
  filter(between(Year, min(y_range), max(y_range))) %>% 
  select(-OpenInterval) %>% 
  mutate(yg = cut(Year, breaks = y_breaks, include.lowest = T),
         yg = as.character(yg),
         Agegr = cut(Age, age_breaks_fert, age_labels_fert, 
                     include.lowest = TRUE, right = F),
         Year = as.numeric(str_sub(yg,2,5)),) %>% 
  group_by(Year, Agegr) %>% 
  summarise(ASFR = mean(ASFR)) %>% 
  ungroup %>% 
  mutate(Source = "HFD")

# Wrangling Socsim data
SocsimF0 <- asfr_wod %>% 
  rename(Agegr = age, ASFR = socsim) %>% 
  mutate(Year = as.numeric(str_sub(year,2,5)),
         Source = "Socsim") %>% 
  select(-year)

## Plotting ASFR from HFD vs Socsim   
bind_rows(HFD0, SocsimF0) %>% 
  filter(!is.na(ASFR)) %>% 
  # filter(Year == 1970) %>% 
  ggplot(aes(x = Agegr, y = ASFR, group = interaction(Year, Source)))+
  geom_line(aes(colour = Source, linetype = Source), size = 1)+
  scale_color_manual(values=c("#30734A", "#CA650D"))+
  theme_bw()


# Calculate TFR from HFD
HFD1 <- HFD %>% 
  filter(between(Year, min(y_range), max(y_range))) %>% 
  group_by(Year) %>% 
  summarise(TFR = sum(ASFR, na.rm=T)) %>% 
  ungroup %>% 
  mutate(Source = "HFD") 

HFD1 %>%
  ggplot(aes(x = Year, y = TFR)) +
  geom_line(size=1, colour = "#30734A") +
  theme_bw()

# Calculate TFR from Socsim
SocsimF <-  asfr_wod %>% 
  group_by(year) %>% 
  summarise(FR = sum(socsim, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(Year = as.numeric(str_sub(year,2,5)), 
         TFR = FR*5,
         Source = "Socsim") 

## Plotting TFR from HFD vs Socsim 
bind_rows(HFD1, SocsimF) %>% 
  ggplot(aes(x = Year, y = TFR, colour = Source)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#30734A", "#CA650D"))+
  theme_bw() 


##  Mortality rates ----

# HMD credentials
# HMD_username <- "Type_here_HMD_username"
# HMD_password <- "Type_here_HMD_password"

# Get female life tables from HMD
ltf <- readHMDweb(CNTRY = "SWE",
                  item = "fltper_5x5",
                  username = HMD_username,
                  password = HMD_password)

# Get male life tables from HMD
ltm <- readHMDweb(CNTRY = "SWE",
                  item = "mltper_5x5",
                  username = HMD_username,
                  password = HMD_password)


# Wrangling HMD life tables 
## We will compare central death rates from life tables (not qx as in the input rates)
HMD <- ltf %>%
  select(Year, Age, mx) %>% 
  mutate(Sex = "Female") %>% 
  bind_rows(ltm %>% 
              select(Year, Age, mx) %>%  
              mutate(Sex = "Male")) %>% 
  mutate(Sex = factor(Sex, levels = c("Male", "Female")),
         Source = "HMD")

## Plotting Central Death Rates vs Socsim 
HMD %>% 
  ggplot(aes(x = Age, y = mx, group = Year))+
  facet_wrap(~Sex) +
  geom_line(aes(colour = Year))+
  scale_y_log10() +
  theme_bw() + 
  scale_color_viridis(option = "A") +
  guides(colour = guide_colourbar(reverse = TRUE)) +
  labs(title = "Central death rates by sex in Sweden (from HMD life tables)", 
       x = "Age", y = "mx") 

# Wrangling Socsim data
SocsimM <- asmr_wod %>% 
  rename(Age = age, mx = socsim) %>% 
  mutate(Year = as.numeric(str_sub(year,2,5)),
         Sex = ifelse(sex == "male", "Male", "Female"),
         Source = "Socsim") %>% 
  select(Year, Age, Sex, mx, Source)

## Plotting ASMR from HMD vs Socsim 
HMD %>%
  filter(between(Year, min(y_range), max(y_range)) & Age <= 95) %>% 
  bind_rows(SocsimM) %>% 
  filter(!is.na(mx) & mx != 0) %>% 
  ggplot(aes(x = Age, y = mx, group = interaction(Year, Source)))+
  facet_wrap(~Sex) +
  geom_line(aes(colour = Source, linetype = Source), size = 1)+
  scale_y_log10() +
  scale_color_manual(values=c("#0C0B7F", "#CA650D"))+
  theme_bw()

#------------------------------------------------------------------------------------------------------
## Trace ancestors of birth cohorts (generations) ----
# Ancestors of a sub-set of birth cohorts born every 50 years from 1800-2000

# Select years to analyse
# periods <- seq(1800, 2000, by = 50)
# 
# # Collect pids of people born at each time point of the period
# bchs <- opop %>% 
#   mutate(Generation = asYr(dob)) %>% 
#   filter(Generation %in% periods) %>% 
#   pull(pid)
# 
# # Get a sample of 10% of people born in the above-mentioned decades
# sample_size <- round(length(bchs)/10)
# bchs_samp <-  sample(bchs, sample_size, replace = F)
# 
# ## Get a glimpse on egos' information
# opop %>% filter(pid %in% bchs_samp)

## Map the function to get the ancestors of multiple egos in a single dfr
# start <- Sys.time()
# ancestors_bchs <- map_dfr(bchs_samp, get_ancestors) 
# end <- Sys.time()
# print(end-start)

# Distribution of egos by birth cohort (Generation)
# ancestors_bchs %>% 
#   filter(pid == ego_id) %>% 
#   mutate(Generation = asYr(dob)) %>% 
#   pull(Generation) %>% 
#   freq()

# Check number of duplicates (summing True values)
# ancestors_bchs %>% pull(pid) %>% duplicated() %>% sum()

#------------------------------------------------------------------------------------------------------