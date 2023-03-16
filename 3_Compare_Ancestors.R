#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Genealogies - Trace and compare subset of ancestors of a given ego(s) 
# U:/SOCSIM/SOCSIM_Genealogies/3_Compare_Ancestors.R

## Trace ancestors of a given ego(s) from a SOCSIM microsimulation for Sweden (1751-2021)
# and compare demographic measures from the whole simulation and the genealogical subsets

# Created by Liliana Calderon on 23-09-2022
# Last modified by Liliana Calderon on 16-03-2023

## NB: To run this code, it is necessary to have already run the script 1_Run_Simulations.R

#------------------------------------------------------------------------------------------------------
## General settings ----
## We use only one of the 10 simulations. 
options(scipen=999999)

## Load packages 
library(tidyverse)
library(ggh4x)  # To facet scales-
library(questionr)

## Read the output .opop file ----

# Load read_opop() function to read the .opop file (written by Diego Alburez-Gutierrez)
# Once this is integrated into the rsocsim package the function might be just called
source("read_opop.R")

# Name of the supervisory file used for the simulation (if not set in the GlobalEnv)
supfile <- "socsim_SWE_marr.sup"

## Randomly choose the simulation seed to use 
load("sims_seeds.rda")
seed <- "13486"
# seed <-  sample(sims_seeds, 1, replace = F) 

# Path of the folder where the simulation results have been saved. 
path_opop <- paste0("sim_results_s",supfile,"_",seed,"/result.opop")

# Read SOCSIM output, using read_opop function. 
opop <- read_opop(path_opop)

#------------------------------------------------------------------------------------------------------
## Trace ancestors of people alive in 2022 as a proxy of current genealogists 

## Load functions to get ancestors
source("Functions_Ancestors.R")

# Pids of people alive at the end of the simulation, i.e. dod == 0, 
# who are older than 18 years old on 01-01-2022, i,e. dob 1913-2003 
egos2022 <- opop %>% 
  mutate(Generation = asYr(dob, last_month, final_sim_year)) %>% 
  filter(dod == 0 & Generation <= final_sim_year-18) %>% 
  pull(pid)

# Get a sample of 1% of people alive in 2022. 
# To be increased after
sample_size <- round(length(egos2022)/100)
egos2022_samp <-  sample(egos2022, sample_size, replace = F)
save(egos2022_samp, file = "egos2022_samp.RData")

## Map the function to get the ancestors of a sample of individuals alive in 2022 (older than 18 years)
start <- Sys.time()
ancestors_egos2022 <- map_dfr(egos2022_samp, get_ancestors)
end <- Sys.time()
print(end-start) # Time difference of 7.810757 hours

# Save the data frame
save(ancestors_egos2022, file = "ancestors_egos2022.RData")

#----------------------------------------------------------------------------------------------------
## Recover age-specific fertility and mortality rates  -----
# Retrieve and compare the rates derived from the whole simulation with those from the genealogical subset,
# of direct ancestors both will and without duplicates

# Prevent scientific notation (useful for the rate calculation)
options(scipen=999999)

## Load functions to recover input age-specific fertility and mortality rates 
source("Functions_Retrieve_Rates.R")

## Calculate ASFR and ASMR for the Whole simulation (seed "13486")

# Retrieve age-specific fertility rates for the whole single simulation 
asfr_whole <- get_asfr_socsim(df = opop,
                              final_sim_year = 2021 , #[Jan-Dec]
                              year_min = 1750, # Closed [
                              year_max = 2020, # Open )
                              year_group = 5, 
                              age_min_fert = 10, # Closed [
                              age_max_fert = 55, # Open )
                              age_group = 5) #[,)
save(asfr_whole, file = "asfr_whole.RData")

# Retrieve age-specific mortality rates for the whole single simulation
asmr_whole <- get_asmr_socsim(df = opop, 
                              final_sim_year = 2021, #[Jan-Dec]
                              year_min = 1750, # Closed
                              year_max = 2020, # Open )
                              year_group = 5,
                              age_max_mort = 110, # Open )
                              age_group = 5) #[,)
save(asmr_whole, file = "asmr_whole.RData")

# Load the data frame with the ancestors of 1% sample of egos alive in 2022
load("ancestors_egos2022.RData")

#  Calculate ASFR and ASMR for genealogical subset of direct ancestors of population alive in 01-01-2022 with duplicates 

# Copy the vector of direct ancestors with duplicates. 
ancestors_egos2022_wd <- ancestors_egos2022 

# Retrieve age-specific fertility rates for the genealogical subset of direct ancestors with duplicates
asfr_wd <- get_asfr_socsim(df = ancestors_egos2022_wd,
                           final_sim_year = 2021 , #[Jan-Dec]
                           year_min = 1750, # Closed [
                           year_max = 2020, # Open )
                           year_group = 5, 
                           age_min_fert = 10, # Closed [
                           age_max_fert = 55, # Open )
                           age_group = 5) #[,)
save(asfr_wd, file = "asfr_wd.RData")

# Retrieve age-specific mortality rates for the genealogical subset of direct ancestor with duplicates
asmr_wd <- get_asmr_socsim(df = ancestors_egos2022_wd,
                           final_sim_year = 2021, #[Jan-Dec]
                           year_min = 1750, # Closed
                           year_max = 2020, # Open )
                           year_group = 5,
                           age_max_mort = 110, # Open )
                           age_group = 5) #[,)
save(asmr_wd, file = "asmr_wd.RData")


## Calculate ASFR and ASMR for genealogical subset of direct ancestors without duplicates 

# Ancestors without duplicates for sample 
ancestors_egos2022_wod <- ancestors_egos2022 %>% distinct(pid, .keep_all = TRUE)

# Retrieve age-specific fertility rates for the genealogical subset of direct ancestors without duplicates
asfr_wod <- get_asfr_socsim(df = ancestors_egos2022_wod,
                            final_sim_year = 2021 , #[Jan-Dec]
                            year_min = 1750, # Closed [
                            year_max = 2020, # Open )
                            year_group = 5, 
                            age_min_fert = 10, # Closed [
                            age_max_fert = 55, # Open )
                            age_group = 5) #[,)
save(asfr_wod, file = "asfr_wod.RData")

# Retrieve age-specific mortality rates for the genealogical subset of direct ancestors without duplicates
asmr_wod <- get_asmr_socsim(df = ancestors_egos2022_wod,
                            final_sim_year = 2021, #[Jan-Dec]
                            year_min = 1750, # Closed
                            year_max = 2020, # Open )
                            year_group = 5,
                            age_max_mort = 110, # Open )
                            age_group = 5) #[,)
save(asmr_wod, file = "asmr_wod.RData")

#----------------------------------------------------------------------------------------------------
## Plot results for the genealogical subsets of direct ancestors with and without duplicates ----
library(viridis)
library(svglite) # To save svg files

## Theme for the graphs
source("Functions_Graphs.R")

# Load ASFR and ASMR for the genealogical subset of direct ancestors with duplicates
load("asfr_wd.RData")
load("asmr_wd.RData")

# Load ASFR and ASMR for the genealogical subset of direct ancestors without duplicates
load("asfr_wod.RData")
load("asmr_wod.RData")

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
load("asfr_whole.RData")
load("asmr_whole.RData")

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
# labs(title = "Age-Specific Fertility Rates in Sweden (1751-2021), 
# retrieved from a SOCSIM simulation and subsets of direct ancestors") 
ggsave(file="Graphs/Socsim_ances_ASFR.jpeg", width=17, height=9, dpi=400)


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
# yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

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
  #labs(title = "Age-Specific Mortality Rates in Sweden (1751-2021), retrieved from a SOCSIM simulation and subset of direct ancestors") 
ggsave(file="Graphs/socsim_ances_ASMR.jpeg", width=17, height=9, dpi=400)


## Final plot combining ASFR and ASMR ----

# Same years to plot than above (in intervals). Change if necessary
# yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

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
  mutate(age = factor(as.character(age), levels = age_levels)) %>% 
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset)))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(aes(colour = Year, linetype = Dataset), linewidth = 1.3)+
  scale_color_manual(values = c("#FC8961", "#B72779", "#51127C"))+
  scale_linetype_manual(values = c("11", "22", "solid")) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(x = "Age")
ggsave(file="Graphs/Final_Socsim_Ances_ASFR_ASMR.jpeg", width=17, height=9, dpi=400)

## Save as .svg file for poster
ggsave(file="Graphs/socsim_ances_ASFR_ASMR.svg", device = "svg", units = "in", width=15, height=8, dpi=400) 

#----------------------------------------------------------------------------------------------------
#### Summary measures: TFR and e0 ----
# Here, we use the rates by single calendar year and 1 year of age

#### Total Fertility Rate ----
# Calculate Total Fertility Rate from asfr 1x1

# Retrieve age-specific fertility rates 1x1 for the whole single simulation 
asfr_whole_1 <- get_asfr_socsim(df = opop,
                                final_sim_year = 2021 , #[Jan-Dec]
                                year_min = 1750, # Closed [
                                year_max = 2020, # Open )
                                year_group = 1, 
                                age_min_fert = 10, # Closed [
                                age_max_fert = 55, # Open )
                                age_group = 1) #[,)
save(asfr_whole_1, file = "asfr_whole_1.RData")

# Retrieve age-specific fertility rates 1x1 for the genealogical subset of direct ancestors with duplicates
asfr_wd_1 <- get_asfr_socsim(df = ancestors_egos2022_wd,
                              final_sim_year = 2021 , #[Jan-Dec]
                              year_min = 1750, # Closed [
                              year_max = 2020, # Open )
                              year_group = 1, 
                              age_min_fert = 10, # Closed [
                              age_max_fert = 55, # Open )
                              age_group = 1) #[,)
save(asfr_wd_1, file = "asfr_wd_1.RData")

# Retrieve age-specific fertility rates 1x1 for the genealogical subset of direct ancestors without duplicates
asfr_wod_1 <- get_asfr_socsim(df = ancestors_egos2022_wod,
                            final_sim_year = 2021 , #[Jan-Dec]
                            year_min = 1750, # Closed [
                            year_max = 2020, # Open )
                            year_group = 1, 
                            age_min_fert = 10, # Closed [
                            age_max_fert = 55, # Open )
                            age_group = 1) #[,)
save(asfr_wod_1, file = "asfr_wod_1.RData")



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

## Plotting TFR from whole SOCSIM simulation and the genealogical subsets of direct ancestors with duplicates
bind_rows(TFR_whole, TFR_wd, TFR_wod) %>% 
  filter(Year >= 1751) %>%
  ggplot(aes(x = Year, y = TFR)) +
  geom_line(aes(colour = Dataset, linetype = Dataset), linewidth = 1.3)+ 
  scale_color_manual(values = c("#FC8961", "#B73779", "#1A7761"))+
  scale_linetype_manual(values = c("11", "22", "solid")) +
  theme_graphs() 
# labs(title = "Total Fertility Rate in Sweden (1751-2021), retrieved from a SOCSIM simulation and genealogical subset of direct ancestors") 
ggsave(file="Graphs/socsim_ances_TFR.jpeg", width=17, height=9, dpi=400)

## Life Expectancy at birth ----
# Calculate life expectancy at birth from asmr 1x1

# Retrieve age-specific mortality rates for the whole simulation
asmr_whole_1 <- get_asmr_socsim(df = opop, 
                              final_sim_year = 2021, #[Jan-Dec]
                              year_min = 1750, # Closed
                              year_max = 2020, # Open )
                              year_group = 1,
                              age_max_mort = 110, # Open )
                              age_group = 1) #[,)
save(asmr_whole_1, file = "asmr_whole_1.RData")


# Retrieve age-specific mortality rates for the genealogical subset of direct ancestors with duplicates
asmr_wd_1 <- get_asmr_socsim(df = ancestors_egos2022_wd,
                             final_sim_year = 2021, #[Jan-Dec]
                             year_min = 1750, # Closed
                             year_max = 2020, # Open )
                             year_group = 1,
                             age_max_mort = 110, # Open )
                             age_group = 1) #[,)
save(asmr_wd_1, file = "asmr_wd_1.RData")

# Retrieve age-specific mortality rates for the genealogical subset of direct ancestors without duplicates
asmr_wod_1 <- get_asmr_socsim(df = ancestors_egos2022_wod,
                              final_sim_year = 2021, #[Jan-Dec]
                              year_min = 1750, # Closed
                              year_max = 2020, # Open )
                              year_group = 1,
                              age_max_mort = 110, # Open )
                              age_group = 1) #[,)
save(asmr_wod_1, file = "asmr_wod_1.RData")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with the asmr df calculated with the get_asmr_socsim()
source("Functions_Life_Table.R")

# Compute life table from asmr 1x1 for Whole SOCSIM simulation
lt_whole <- lt_socsim(asmr_whole_1)
save(lt_whole, file = "lt_whole.RData")

# Compute life table from asmr 1x1 for the genealogical subset of direct ancestors with duplicates
lt_wd <- lt_socsim(asmr_wd_1)
save(lt_wd, file = "lt_wd.RData")

# Compute life table from asmr 1x1 for the genealogical subset of direct ancestors without duplicates
lt_wod <- lt_socsim(asmr_wod_1)
save(lt_wod, file = "lt_wod.RData")

#CHECK
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
ggsave(file="Graphs/socsim_ances_e0.jpeg", width=17, height=9, dpi=400)


#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

## Plotting TFR and e0 (for females) from HFD/HMD vs SOCSIM 

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
  mutate(Rate = factor(Rate, levels = c("TFR", "e0"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(aes(color = Dataset, linetype = Dataset), linewidth = 1.2) +
  scale_color_manual(values = c("#771A30", "#331A77", "#1A7761"))+
  scale_linetype_manual(values = c("11", "22", "solid")) +
  theme_graphs()
# labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2020), retrieved from HFD, HMD and 10 SOCSIM simulation outputs") + 

# Save the plot
ggsave(file="Graphs/Final_Socsim_Ances_TFR_e0.jpeg", width=17, height=9, dpi=400)