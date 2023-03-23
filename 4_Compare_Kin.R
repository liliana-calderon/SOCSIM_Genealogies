#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Sweden - Trace and compare a subset of kin of a given ego(s) 
# U:/SOCSIM/SOCSIM_Genealogies/4_Compare_Kin.R

## Trace relevant ascending and lateral kin up to the 4th generation
# of a given ego(s) from a SOCSIM microsimulation for Sweden (1751-2021)
# and compare demographic measures from the whole simulation and the subsets of family trees

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

## Load functions to get ascending and lateral kin up to the 4th degree of consanguinity
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

# Load same opop sample used to trace the ancestors in 3_Compare_Ancestors.R
# This is a 10% sample of people alive at the end of the simulation, i.e. dod == 0, 
# who are older than 18 years old on 01-01-2022, i,e. dob 1913-2003 

load("egos2022_samp_10.RData")

## Map the function to get relevant kin of a sample of individuals alive in 2022 (older than 18 years)
start <- Sys.time()
kin_egos2022_10 <- map_dfr(egos2022_samp, get_kin) %>%
  left_join(select(opop, c(pid, fem, dob, dod, mom, marid, mstat)), by = "pid")
end <- Sys.time()
print(end-start) # Time difference of 11.66745 hours

# Save the data frame
save(kin_egos2022_10, file = "kin_egos2022_10.RData")

#----------------------------------------------------------------------------------------------------
## Recover age-specific fertility and mortality rates  -----
# Retrieve and compare rates derived from the whole simulation with those from the subset of family trees
# of ascending and lateral kin up to the 4th generation

load("kin_egos2022_10.RData")

##  Calculate ASFR and ASMR for the subset of "extended" family trees up to the 4th degree of consanguinity:  

# Ascending and lateral kin up to the 4th degree for sample of egos alive in 2022 without duplicates
kin_egos2022_ext <- kin_egos2022_10 %>% distinct(pid, .keep_all = TRUE)

# Check minimum year of birth to define year_min
kin_egos2022_ext %>% 
  mutate(last_month = max(dob),
         final_sim_year = 2021,
         Birth_Year = asYr(dob, last_month, final_sim_year)) %>% 
  pull(Birth_Year) %>% 
  min()

# Since the minimum year of birth here is 1752, 
# the minimum year_min for the fertility rates calculation should be 1762 (rounded to 1765) 
# as age_min_fert = 10 women could be counted in the first age group after age 10. 
# For mortality, 1755 could be the year_min

# Retrieve age-specific fertility rates for the subset of "extended" family trees without duplicates
asfr_ext <- get_asfr_socsim(df = kin_egos2022_ext,
                           final_sim_year = 2021 , #[Jan-Dec]
                           year_min = 1765, # Closed [
                           year_max = 2020, # Open )
                           year_group = 5, 
                           age_min_fert = 10, # Closed [
                           age_max_fert = 55, # Open )
                           age_group = 5) #[,)

save(asfr_ext, file = "asfr_ext.RData")

# Retrieve age-specific mortality rates for the subset of "extended" family trees without duplicates
asmr_ext <- get_asmr_socsim(df = kin_egos2022_ext,
                           final_sim_year = 2021, #[Jan-Dec]
                           year_min = 1755, # Closed
                           year_max = 2020, # Open )
                           year_group = 5,
                           age_max_mort = 110, # Open )
                           age_group = 5) #[,)
save(asmr_ext, file = "asmr_ext.RData")


##  Calculate ASFR and ASMR for the subset of "direct" family trees up to the 4th degree of consanguinity

# Only direct ascending kin up to the 4th degree for sample of egos alive in 2022, without duplicates
# This is similar to direct ancestors up to the 5th generation. 
ancestors_5 <- c("ego", "m", "f", "mm","mf", "fm", "ff", 
                 "mmm", "mmf", "mfm", "mff", "fmm", "fmf", "ffm", "fff", 
                 "mmmm", "mmmf", "mmfm", "mmff", "mfmm", "mfmf", "mffm", "mfff", 
                 "fmmm", "fmmf", "fmfm", "fmff", "ffmm", "ffmf", "fffm", "ffff")
kin_egos2022_dir <- kin_egos2022_ext %>% 
  filter(kin_type %in% ancestors_5)

# Check minimum year of birth to define year_min
kin_egos2022_dir %>% 
  mutate(last_month = max(dob),
         final_sim_year = 2021, 
         Birth_Year = asYr(dob, last_month, final_sim_year)) %>% 
  pull(Birth_Year) %>% 
  min() # 1770

# The minimum year_min for the fertility rates calculation should be 1785
# as age_min_fert = 10 women could be counted in the first age group after age 10. 

# Retrieve age-specific fertility rates for the subset of "direct" family trees without duplicates
asfr_dir <- get_asfr_socsim(df = kin_egos2022_dir,
                            final_sim_year = 2021 , #[Jan-Dec]
                            year_min = 1785, # Closed [
                            year_max = 2020, # Open )
                            year_group = 5, 
                            age_min_fert = 10, # Closed [
                            age_max_fert = 55, # Open )
                            age_group = 5) #[,)
save(asfr_dir, file = "asfr_dir.RData")

# Retrieve age-specific mortality rates for the subset of "direct" family trees without duplicates
asmr_dir <- get_asmr_socsim(df = kin_egos2022_dir,
                            final_sim_year = 2021, #[Jan-Dec]
                            year_min = 1770, # Closed
                            year_max = 2020, # Open )
                            year_group = 5,
                            age_max_mort = 110, # Open )
                            age_group = 5) #[,)
save(asmr_dir, file = "asmr_dir.RData")

#----------------------------------------------------------------------------------------------------
## Plot results for the genealogical subsets of ascending and lateral kin up to the 4th generation ----

# Load ASFR and ASMR for the subset of "direct" family trees without duplicates
load("asfr_dir.RData")
load("asmr_dir.RData")

# Load ASFR and ASMR for the subset of "extended" family trees without duplicates
load("asfr_ext.RData")
load("asmr_ext.RData")

# Choose years to plot (in intervals).
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_ext$age)


## ASFR and ASMR (for women) for the subset of "extended" family trees without duplicates
bind_rows(asfr_ext %>% 
            mutate(rate = "ASFR",
                   sex = "female"),
          asmr_ext %>% 
            mutate(rate = "ASMR") %>% 
            filter(sex == "female")) %>% 
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(socsim !=0 & !is.infinite(socsim) & !is.nan(socsim)) %>% 
  filter(year %in% yrs_plot) %>% 
  mutate(age = factor(as.character(age), levels = age_levels)) %>% 
  ggplot(aes(x = age, y = socsim, group = year, colour = year)) +
  geom_line(linewidth = 1) +
  facet_wrap(. ~ rate, scales = "free") + 
  theme_graphs()  +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10"))) +
  scale_color_viridis(option = "F", discrete = T, direction = -1)+
  labs(x = "Age", y = "Estimate")

## ASFR and ASMR (for women) for the subset of "direct" family trees without duplicates
bind_rows(asfr_dir %>% 
            mutate(rate = "ASFR",
                   sex = "female"),
          asmr_dir %>% 
            mutate(rate = "ASMR") %>% 
            filter(sex == "female")) %>% 
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(socsim !=0 & !is.infinite(socsim) & !is.nan(socsim)) %>% 
  filter(year %in% yrs_plot) %>% 
  mutate(age = factor(as.character(age), levels = age_levels)) %>% 
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
## Comparison of a whole SOCSIM simulation with family trees subsets of ascending and lateral kin ----

# Load ASFR and ASMR for the whole single simulation (seed "13486"), used in 3_Compare_Ancestors.R
load("asfr_whole.RData")
load("asmr_whole.RData")

#### Age-specific Fertility Rates ----

# Whole SOCSIM simulation
asfr_whole2 <- asfr_whole %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Whole_simulation", 
         Rate = "ASFR") 

# "Extended" family trees without duplicates
asfr_ext2 <- asfr_ext %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Trees_extended", 
         Rate = "ASFR") 

# "Direct" family trees without duplicates
asfr_dir2 <- asfr_dir %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Trees_direct",
         Rate = "ASFR") 

## Plot ASFR from whole SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2, asfr_ext2, asfr_dir2) %>% 
  filter(year %in% yrs_plot & !is.nan(ASFR)) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset)))+
  geom_line(aes(colour = year, linetype = Dataset), linewidth = 1.2)+ 
  scale_color_viridis(option = "D", discrete = T, direction = -1) +
  scale_linetype_manual(values = c("11", "22", "solid")) +
  theme_graphs()
# labs(title = "Age-Specific Fertility Rates in Sweden (1751-2021), 
# retrieved from a SOCSIM simulation and subsets of "extended" and direct" family trees") 
ggsave(file="Graphs/Socsim_Trees_ASFR.jpeg", width=17, height=9, dpi=400)


## Age-Specific Mortality rates ----

# Whole SOCSIM simulation
asmr_whole2 <- asmr_whole %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Whole_simulation",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# "Extended" family trees without duplicates
asmr_ext2 <- asmr_ext %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Trees_extended",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)

# "Direct" family trees without duplicates
asmr_dir2 <- asmr_dir %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Trees_direct",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)


## Plotting ASMR from whole SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates

# Same years to plot than above (in intervals). Change if necessary
# yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asmr_whole2, asmr_ext2, asmr_dir2) %>% 
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
#labs(title = "Age-Specific Mortality Rates in Sweden (1751-2021), 
# retrieved from a SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates") 
ggsave(file="Graphs/Socsim_Trees_ASMR.jpeg", width=17, height=9, dpi=400)


## Final plot combining ASFR and ASMR ----

# Same years to plot than above (in intervals). Change if necessary
# yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_whole2$age)

## Ploting ASFR and ASMR (for females) from whole SOCSIM simulation and subsets of "extended" and direct" family trees
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_ext2 %>% rename(Estimate = ASFR),
          asfr_dir2 %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx),
            asmr_ext2 %>% rename(Estimate = mx),
            asmr_dir2 %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels)) %>% 
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset)))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(aes(colour = Year, linetype = Dataset), linewidth = 1.3)+
  scale_color_manual(values = c("#FC8961", "#B72779", "#2779B7"))+
  scale_linetype_manual(values = c("11", "22", "solid")) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(x = "Age")
ggsave(file="Graphs/Final_Socsim_Trees_ASFR_ASMR.jpeg", width=17, height=9, dpi=400)

## Save as .svg file for poster
# ggsave(file="Graphs/Socsim_Trees_ASFR_ASMR.svg", device = "svg", units = "in", width=15, height=8, dpi=400) 

#----------------------------------------------------------------------------------------------------
#### Summary measures: TFR and e0 ----
# Here, we use the rates by single calendar year and 1 year of age

#### Total Fertility Rate ----
# Calculate Total Fertility Rate from asfr 1x1

# Load ASFR 1x1 for the whole single simulation (seed "13486"), used in 3_Compare_Ancestors.R
load("asfr_whole_1.RData")

# Retrieve age-specific fertility rates 1x1 for the subset of "extended" family trees without duplicates
asfr_ext_1 <- get_asfr_socsim(df = kin_egos2022_ext,
                              final_sim_year = 2021 , #[Jan-Dec]
                              year_min = 1765, # Closed [
                              year_max = 2020, # Open )
                              year_group = 1, 
                              age_min_fert = 10, # Closed [
                              age_max_fert = 55, # Open )
                              age_group = 1) #[,)
save(asfr_ext_1, file = "asfr_ext_1.RData")
load("asfr_ext_1.RData")

# Retrieve age-specific fertility rates 1x1 for the subset of "direct" family trees without duplicates
asfr_dir_1 <- get_asfr_socsim(df = kin_egos2022_dir,
                              final_sim_year = 2021 , #[Jan-Dec]
                              year_min = 1785, # Closed [
                              year_max = 2020, # Open )
                              year_group = 1, 
                              age_min_fert = 10, # Closed [
                              age_max_fert = 55, # Open )
                              age_group = 1) #[,)
save(asfr_dir_1, file = "asfr_dir_1.RData")
load("asfr_dir_1.RData")

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

# "Extended" family trees subset without duplicates
TFR_ext <- asfr_ext_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "Trees_extended",
         Rate = "TFR", 
         sex = "female") 

# "Direct" family trees without duplicates
TFR_dir <- asfr_dir_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "Trees_direct",
         Rate = "TFR", 
         sex = "female")

## Plott TFR from whole SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates
bind_rows(TFR_whole, TFR_ext, TFR_dir) %>% 
  filter(Year >= 1751) %>%
  ggplot(aes(x = Year, y = TFR)) +
  geom_line(aes(colour = Dataset, linetype = Dataset), linewidth = 1.3)+ 
  scale_color_manual(values = c("#5B4792", "#621528", "#1A7761"))+
  scale_linetype_manual(values = c("22", "11", "solid")) +
  theme_graphs() 
# labs(title = "Total Fertility Rate in Sweden (1751-2021), 
#retrieved from a SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates") 
ggsave(file="Graphs/Socsim_Trees_TFR.jpeg", width=17, height=9, dpi=400)

## Life Expectancy at birth ----
# Calculate life expectancy at birth from asmr 1x1

# Load age-specific mortality rates 1x1 for the whole single simulation (seed "13486"),
load("asmr_whole_1.RData")
# Load life table from asmr 1x1 for the whole single simulation (seed "13486"), used in 3_Compare_Ancestors.R
load("lt_whole.RData")

# Retrieve age-specific mortality rates for the subset of "extended" family trees without duplicates
asmr_ext_1 <- get_asmr_socsim(df = kin_egos2022_ext,
                              final_sim_year = 2021, #[Jan-Dec]
                              year_min = 1755, # Closed
                              year_max = 2020, # Open )
                              year_group = 1,
                              age_max_mort = 110, # Open )
                              age_group = 1) #[,)
save(asmr_ext_1, file = "asmr_ext_1.RData")

# Retrieve age-specific mortality rates for the subset of "direct" family trees without duplicates
asmr_dir_1 <- get_asmr_socsim(df = kin_egos2022_dir,
                              final_sim_year = 2021, #[Jan-Dec]
                              year_min = 1770, # Closed
                              year_max = 2020, # Open )
                              year_group = 1,
                              age_max_mort = 110, # Open )
                              age_group = 1) #[,)
save(asmr_dir_1, file = "asmr_dir_1.RData")

## In the first years of the subset of family trees, at almost all ages 
# rates have values of 0, Infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
# Hence, the life table cannot be computed and must be filter after some year
# where rates are available for the whole age range. 

# Load asmr 1x1 for the subset of "extended" family trees. 
load("asmr_ext_1.RData")

# Filter after year without warnings in filter(between(rn, 1, max(which(mx > 0))))
asmr_ext_1_filt <- asmr_ext_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  filter(Year > 1808)

# Compute life table from asmr 1x1 for the subset of "extended" family trees without duplicates
lt_ext <- lt_socsim(asmr_ext_1_filt)
save(lt_ext, file = "lt_ext.RData")

# Load asmr 1x1 for the subset of "direct" family trees
load("asmr_dir_1.RData")

# Filter after year without warnings in filter(between(rn, 1, max(which(mx > 0))))
asmr_dir_1_filt <- asmr_dir_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  filter(Year > 1826)

# Compute life table from asmr 1x1 for the subset of "direct" family trees without duplicates
lt_dir <- lt_socsim(asmr_dir_1_filt)
save(lt_dir, file = "lt_dir.RData")


# Wrangle data for plotting

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

# "Extended" family trees subset without duplicates
lt_ext2 <- lt_ext %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Trees_extended",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

# "Direct" family trees without duplicates
lt_dir2 <- lt_dir %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Trees_direct",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

bind_rows(lt_whole2, lt_ext2, lt_dir2) %>% 
  filter(Age == 0 & Year %in% year_range_mort_1) %>%
  mutate(Sex = ifelse(sex == "female", "Female", "Male")) %>% 
  ggplot(aes(x = Year, y = ex, group = Dataset))+
  geom_line(aes(colour = Dataset, linetype = Dataset), linewidth = 1.3)+
  scale_color_manual(values = c("#5B4792", "#621528", "#1A7761"))+
  scale_linetype_manual(values = c("22", "11", "solid")) +
  facet_wrap(~Sex) +
  theme_graphs() +
  labs(y = "e0") 
 #  labs(title = "Life expectancy at birth in Sweden (e0), 1751-2020, retrieved from a SOCSIM simulation 
       # subsets of "extended" and direct" family trees without duplicates")
ggsave(file="Graphs/socsim_Trees_e0.jpeg", width=17, height=9, dpi=400)


#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

## TFR and e0 (for females) from whole SOCSIM simulation and subsets of "extended" and direct" family trees 

bind_rows(TFR_whole %>% 
            rename(Estimate = TFR), 
          TFR_ext %>% 
            rename(Estimate = TFR), 
          TFR_dir %>%           
            rename(Estimate = TFR)) %>%  
  bind_rows(lt_whole2 %>% 
              rename(Estimate = ex) %>% 
              filter(Age == 0),
            lt_ext2 %>% 
              rename(Estimate = ex) %>% 
              filter(Age == 0),
            lt_dir2 %>% 
              rename(Estimate = ex) %>% 
              filter(Age == 0)) %>% 
  filter(sex == "female") %>% 
  mutate(Rate = factor(Rate, levels = c("TFR", "e0"))) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(aes(color = Dataset, linetype = Dataset), linewidth = 1.2) +
  scale_color_manual(values = c("#5B4792", "#621528", "#1A7761"))+
  scale_linetype_manual(values = c("22", "11", "solid")) +
  theme_graphs()
# labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2020), retrieved 
# from SOCSIM microsimulation and subsets of "extended" and direct" family trees ") + 

# Save the plot
ggsave(file="Graphs/Final_Socsim_Trees_TFR_e0.jpeg", width=17, height=9, dpi=400)

#----------------------------------------------------------------------------------------------------
## Sex Ratio at Birth and Infant Mortality Rate
# The Functions_Retrieve_Rates.R must be called to use the asYr() function

## Define years of not set in the Global Environment
final_sim_year <- 2021 #[Jan-Dec]
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

# Sex Ratio at Birth by year for "Extended" family trees subset without duplicates
SRB_ext <- ancestors_egos2022_ext %>% 
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
         Dataset = "Trees_extended") 

# Sex Ratio at Birth by year for subset of "direct" family trees without duplicates
SRB_dir <- ancestors_egos2022_dir %>% 
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
         Dataset = "Trees_direct") 

# Plotting SRB
bind_rows(SRB_whole, SRB_ext, SRB_dir) %>% 
  filter(Year >= 1751 & !is.na(SRB)) %>% # No births after 2003
  ggplot(aes(x = Year, y = SRB, group = Dataset, color = Dataset, linetype = Dataset))+
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#771A30", "#331A77", "#1A7761"))+ 
  scale_linetype_manual(values = c("11", "22", "solid")) +
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

# Births by year from "Extended" family trees subset without duplicates
Births_ext <- ancestors_egos2022_ext %>% 
  mutate(Year = asYr(dob, last_month, final_sim_year)) %>% 
  filter(Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>% 
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Trees_extended", 
         Event = "Births")

# Births by year from subset of "direct" family trees without duplicates
Births_dir <- ancestors_egos2022_dir %>% 
  mutate(Year = asYr(dob, last_month, final_sim_year)) %>% 
  filter(Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>% 
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Trees_direct", 
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

# Deaths below age 1 (0-11 months) from "Extended" family trees subset without duplicates
Deaths_0_ext <- ancestors_egos2022_ext %>% 
  filter(dod != 0) %>% 
  mutate(age_death_months = dod-dob,
         Year = asYr(dod, last_month, final_sim_year)) %>% 
  filter(age_death_months < 12 & Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>% 
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Trees_extended", 
         Event = "Deaths")

# Deaths below age 1 (0-11 months) from subset of "direct" family trees without duplicates
# There should be no infant mortality in this subset, but let's double check it
Deaths_0_dir <- ancestors_egos2022_dir %>% 
  filter(dod != 0) %>% 
  mutate(age_death_months = dod-dob,
         Year = asYr(dod, last_month, final_sim_year)) %>% 
  filter(age_death_months < 12 & Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>% 
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Trees_direct", 
         Event = "Deaths")

# Calculate and Plot Infant Mortality Rate (IMR)
IMR <- bind_rows(Births_whole, Births_ext, Births_dir, Deaths_0_whole, Deaths_0_ext, Deaths_0_dir) %>%
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
bind_rows(SRB_whole, SRB_ext, SRB_dir) %>% 
  select(Year, Dataset, SRB) %>% 
  full_join(IMR %>% select(Year, Dataset, IMR), by = c("Year", "Dataset")) %>% 
  pivot_longer(SRB:IMR, names_to = "Measure", values_to = "Value") %>% 
  filter(Year >= 1751 & !is.na(Value)) %>%
  mutate(Measure = factor(Measure, levels = c("SRB", "IMR"))) %>%
  ggplot(aes(x = Year, y = Value, group = Dataset, color = Dataset, linetype = Dataset))+
  facet_wrap(. ~ Measure, scales = "free") + 
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#771A30", "#331A77", "#1A7761"))+
  scale_linetype_manual(values = c("11", "22", "solid")) +
  theme_graphs() +
  facetted_pos_scales(y = list(SRB = scale_y_continuous(limits=c(0.7, 1.3)),
                               IMR =  scale_y_continuous())) +
  theme_graphs() +
  theme(axis.title.y = element_blank())
ggsave(file="Graphs/Socsim_Trees_SRB_IMR.jpeg", width=17, height=9, dpi=400)

#----------------------------------------------------------------------------------------------------
## Births and Deaths by year from whole simulation and subsets of "extended" and direct" family trees  -----

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

# Death counts by year from subset of "Extended" family trees without duplicates
Deaths_ext <- ancestors_egos2022_ext %>% 
  filter(dod != 0) %>% 
  mutate(Year = asYr(dod, last_month, final_sim_year)) %>% 
  filter(Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>%
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Trees_extended", 
         Event = "Deaths")

# Death counts by year from subset of "direct" family trees without duplicates
Deaths_dir <- ancestors_egos2022_dir %>% 
  filter(dod != 0) %>% 
  mutate(Year = asYr(dod, last_month, final_sim_year)) %>% 
  filter(Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>%
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Trees_direct", 
         Event = "Deaths")

# Plotting birth and death counts together. 
bind_rows(Births_whole, Births_ext, Births_dir, Deaths_whole, Deaths_ext, Deaths_dir) %>% 
  filter(Year >= 1751) %>%
  ggplot(aes(x = Year, y = n, group = Dataset, color = Dataset, linetype = Dataset))+
  facet_wrap(. ~ Event) + 
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#771A30", "#331A77", "#1A7761"))+ 
  scale_linetype_manual(values = c("11", "22", "solid")) +
  theme_graphs() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(file="Graphs/Socsim_Trees_Births_Deaths.jpeg", width=17, height=9, dpi=400)
# The absolute values are meaningless