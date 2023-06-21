#------------------------------------------------------------------------------------------------------
# SOCSIM - SOCSIM Sweden - Trace and compare a subset of kin of a given ego(s) 
# U:/SOCSIM/SOCSIM_Genealogies/4_Compare_Kin.R

## Trace relevant ascending and lateral kin up to the 4th generation
# of a given ego(s) from a SOCSIM microsimulation for Sweden (1751-2022)
# and compare demographic measures from the whole simulation and the subsets of family trees

# Created by Liliana Calderon on 13-04-2022
# Last modified by Liliana Calderon on 21-06-2023

## NB: To run this code, it is necessary to have already run the script 1_Run_Simulations.R

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

## Load theme for the graphs
source("Functions_Graphs.R")

# Load function to calculate life table from asmr 1x1
# Currently, it only works with asmr calculated with rsocsim::estimate_mortality_rates()
source("Functions_Life_Table.R")

## Load function to re-code the kin_type for the ancestors
source("Functions_Ancestors.R")

#------------------------------------------------------------------------------------------------------
## Trace kin (up to 4th degree) of people alive in 2023 ----

# Load saved list with opop from 10 simulations, generated in 1_Run_Simulations.R
load("sims_opop.RData")
# Load saved list with omar from 10 simulations, generated in 1_Run_Simulations.R
load("sims_omar.RData")

# Load same sample for each opop used to trace the ancestors in 3_Compare_Ancestors.R
# These are a 10% sample of people alive at the end of the simulation, i.e. dod == 0, 
# who are older than 18 years old on 01-01-2023, i,e. dob 1914-2004 
load("Subsets/egos_samp_10.RData")

## Retrieve the relatives of each simulation sample of egos alive in 2023
kin_10_lc <- pmap_dfr(list(sims_opop, sims_omar, egos_samp_10),
                             ~ rsocsim::retrieve_kin(opop = ..1, omar = ..2, pid = ..3,
                                                     KidsOf = KidsOf,
                                                     extra_kintypes = c("unclesaunts", "niblings", "gunclesaunts", "firstcousins"),
                                                     kin_by_sex = F),
                             .id = "Sim_id")
# Time difference of 56.14612 mins for 10 simulations with initial population of 5000, with hetfert0

# Select relevant variables from each opop and add Sim_id to allow merging with the kin_10 df 
sims_opop_id <- map_dfr(sims_opop, ~ select(.x, c(pid, fem, dob, dod, mom, marid, mstat)), 
                        .id = "Sim_id")

# Format the kin data,  put into long format and add columns from each opop
# Also solve some duplication problems from the function
kin_10 <- kin_10_lc %>% 
  mutate(ego_id = unlist(egos_samp_10)) %>% 
  pivot_longer(-c(Sim_id, ego_id), names_to = "kin_type", values_to = "pid") %>% 
  unnest_longer(kin_type:pid) %>% 
  # The function in the package still have some duplicates that need to be removed
  # With the kin_by_sex = F option, we get nieces and nephews in addition to niblings
  filter(!is.na(pid) & !kin_type %in% c("nieces", "nephews")) %>% 
  # Sometimes, a pid is included as more than one kin type, 
  # e.g. ggparents and gunclesaunts or siblings and firstcousins of the same ego
  ## because parents were included as (gg)parents and (g)unclesaunts
  distinct(Sim_id, ego_id, pid, .keep_all = TRUE) %>% 
  left_join(sims_opop_id, by = c("Sim_id", "pid"))

# Save the data frame with the relatives of 10 simulations samples
save(kin_10, file = "Subsets/kin_10.RData")

#------------------------------------------------------------------------------------------------------
## Add direct ancestors of sample of egos alive in 2023 ----

# Load the data frame with the ancestors of the 10 simulations samples of egos alive in 2023
load("Subsets/ancestors_10.RData")

## Re-code the kin_type for the ancestors and remove duplicates (from the different family trees)
ancestors_10 <- ancestors_10 %>% 
  recode_ancestors() %>%
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all = TRUE) %>% 
  ungroup() 

# Merge the two data frames and remove kin types included in both (i.e.g, parents, gparents, ggparents)
anc_kin_10 <- bind_rows(ancestors_10, kin_10) %>% 
  group_by(Sim_id, ego_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup()

# Save the data frame with the ancestors and relatives of 10 simulations samples
save(anc_kin_10, file = "Subsets/anc_kin_10.RData")

#----------------------------------------------------------------------------------------------------
## Age-Specific Fertility and Mortality rates, 5x5  -----
# Retrieve and compare rates derived from the whole simulation with those from the genealogical subsets
# of ascending and lateral kin up to the 4th generation

load("Subsets/anc_kin_10.RData")

# Keep only unique pids for each simulation and create a list of data frames with opop of ancestors and kin
anc_kin_ext <- anc_kin_10 %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all= TRUE) %>% 
  ungroup() %>%  
  split(.$Sim_id)

##  Estimate ASFR and ASMR for the subset of "extended" family trees  

# Keep only unique pids for each simulation and create a list of data frames containing opop of relatives
kin_ext <- kin_10 %>% 
  group_by(Sim_id) %>% 
  distinct(pid, .keep_all = TRUE) %>% 
  ungroup() %>% 
  split(.$Sim_id)

# Estimate age-specific fertility rates for the subset of direct ancestors and lateral kin without duplicates
asfr_ext <- map_dfr(anc_kin_ext, ~ estimate_fertility_rates(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed [
                                                            year_max = 2020, # Open )
                                                            year_group = 5, 
                                                            age_min_fert = 10, # Closed [
                                                            age_max_fert = 55, # Open )
                                                            age_group = 5), # [,)
                        .id = "Sim_id") 
save(asfr_ext, file = "Measures/asfr_ext.RData")

# Estimate age-specific mortality rates for the subset of direct ancestors and lateral kin without duplicates
asmr_ext <- map_dfr(anc_kin_ext, ~ estimate_mortality_rates(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed
                                                            year_max = 2020, # Open )
                                                            year_group = 5,
                                                            age_max_mort = 110, # Open )
                                                            age_group = 5), # [,)
                    .id = "Sim_id") 
save(asmr_ext, file = "Measures/asmr_ext.RData")

#----------------------------------------------------------------------------------------------------
## Plot estimates from genealogical subsets of ascending and lateral kin----

# Load ASFR and ASMR for the subset of "direct" family trees without duplicates
load("Measures/asfr_dir_wod.RData")
load("Measures/asmr_dir_wod.RData")

# Load ASFR and ASMR for the subset of "extended" family trees without duplicates
load("Measures/asfr_ext.RData")
load("Measures/asmr_ext.RData")

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
bind_rows(asfr_dir_wod %>% 
            mutate(rate = "ASFR",
                   sex = "female"),
          asmr_dir_wod %>% 
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
## Comparison of whole simulation with subsets of ascending and lateral kin ----

# Load ASFR and ASMR for the whole single simulation, used in 3_Compare_Ancestors.R
load("Measures/asfr_whole.RData")
load("Measures/asmr_whole.RData")

# ASFR ----

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
asfr_dir_wod2 <- asfr_dir_wod %>% 
  rename(ASFR = socsim) %>% 
  mutate(Dataset = "Trees_direct",
         Rate = "ASFR") 

## Plot ASFR from whole SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates

# Same years to plot than above (in intervals). Change if necessary
yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asfr_whole2, asfr_ext2, asfr_dir_wod2) %>% 
  filter(year %in% yrs_plot & !is.nan(ASFR)) %>% 
  ggplot(aes(x = age, y = ASFR, group = interaction(year, Dataset)))+
  geom_line(aes(colour = year, linetype = Dataset), linewidth = 1.2)+ 
  scale_color_viridis(option = "D", discrete = T, direction = -1) +
  scale_linetype_manual(values = c("22", "11", "solid")) +
  theme_graphs()
# labs(title = "Age-Specific Fertility Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and subsets of "extended" and direct" family trees") 
ggsave(file="Graphs/Socsim_Exp2_ASFR.jpeg", width=17, height=9, dpi=200)


# ASMR ----

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
asmr_dir_wod2 <- asmr_dir_wod %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Dataset = "Trees_direct",
         Rate = "ASMR") %>% 
  select(year, age, Sex, mx = socsim, Dataset, Rate)


## Plotting ASMR from whole SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates

# Same years to plot than above (in intervals). Change if necessary
 yrs_plot <- c("[1800,1805)", "[1900,1905)", "[2000,2005)") 

bind_rows(asmr_whole2, asmr_ext2, asmr_dir_wod2) %>% 
  rename(Year = year) %>% 
  filter(Year %in% yrs_plot) %>%  
  # Some ages can have rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>% 
  ggplot(aes(x = age, y = mx, group = interaction(Year, Dataset)))+
  facet_wrap(~Sex) +
  geom_line(aes(colour = Year, linetype = Dataset), linewidth = 1.2)+ 
  scale_color_viridis(option = "G", discrete = T, direction = -1) +
  scale_linetype_manual(values = c("22", "11", "solid")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_log10() +
  theme_graphs()
#labs(title = "Age-Specific Mortality Rates in Sweden (1751-2022), 
# retrieved from a SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates") 
ggsave(file="Graphs/Socsim_Exp2_ASMR.jpeg", width=17, height=9, dpi=200)


#----------------------------------------------------------------------------------------------------
## Final plot combining ASFR and ASMR ----

# Years to plot limited to  two years
yrs_plot2 <- c("[1900,1905)", "[2000,2005)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_whole2$age)

## Ploting ASFR and ASMR (for females) from whole SOCSIM simulation and subsets of "extended" and direct" family trees
bind_rows(asfr_whole2 %>% rename(Estimate = ASFR), 
          asfr_ext2 %>% rename(Estimate = ASFR),
          asfr_dir_wod2 %>% rename(Estimate = ASFR)) %>%  
  mutate(Sex = "Female") %>%  
  bind_rows(asmr_whole2 %>% rename(Estimate = mx),
            asmr_ext2 %>% rename(Estimate = mx),
            asmr_dir_wod2 %>% rename(Estimate = mx)) %>% 
  rename(Year = year) %>% 
  filter(Sex == "Female" & Year %in% yrs_plot) %>%
  # There can be rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
  mutate(age = factor(as.character(age), levels = age_levels),
         Rate = ifelse(Rate == "ASFR", "Age-Specific Fertility Rates", 
                       "Age-Specific Mortality Rates"),
         Dataset = case_when(Dataset == "Trees_extended" ~ "Experiment 2 with lateral kin",
                             Dataset == "Trees_direct" ~ "Experiment 2 without lateral kin",
                             TRUE ~ "Whole simulation")) %>% 
  ggplot(aes(x = age, y = Estimate, group = interaction(Year, Dataset), colour = Year))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_line(linewidth = 1.3)+
  geom_point(aes(shape = Dataset), size = 7)+
  scale_color_manual(values = c("#B72779", "#2779B7"))+
  scale_shape_manual(values = c(18, 25, 46)) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_graphs() +
  labs(x = "Age")
ggsave(file="Graphs/Final_Socsim_Exp2_ASFR_ASMR.jpeg", width=17, height=9, dpi=200)

## Save as .svg file for poster
# ggsave(file="Graphs/Socsim_Exp2_ASFR_ASMR.svg", device = "svg", units = "in", width=15, height=8, dpi=200) 

#----------------------------------------------------------------------------------------------------
## Summary measures: TFR and e0 ----
# Here, we use the rates by 1 year age group and 1 calendar year

# Total Fertility Rate ----
# Estimate Total Fertility Rate from asfr 1x1

# Estimate age-specific fertility rates 1x1 for the subset of direct ancestors and lateral kin without duplicates
asfr_ext_1 <- map_dfr(anc_kin_ext, ~ estimate_fertility_rates(opop = .x,
                                                            final_sim_year = 2022, #[Jan-Dec]
                                                            year_min = 1750, # Closed [
                                                            year_max = 2023, # Open )
                                                            year_group = 1, 
                                                            age_min_fert = 10, # Closed [
                                                            age_max_fert = 55, # Open )
                                                            age_group = 1), # [,)
                    .id = "Sim_id") 
save(asfr_ext_1, file = "Measures/asfr_ext_1.RData")

# Load ASFR 1x1 for the whole single simulation (seed "13486"), used in 3_Compare_Ancestors.R
load("Measures/asfr_whole_1.RData")
load("Measures/asfr_ext_1.RData")
load("Measures/asfr_dir_wod_1.RData")

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
TFR_dir_wod <- asfr_dir_wod_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  group_by(Year) %>% 
  summarise(TFR = sum(socsim)*age_group_fert) %>%
  ungroup() %>% 
  mutate(Dataset = "Trees_direct",
         Rate = "TFR", 
         sex = "female")

## Plott TFR from whole SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates
bind_rows(TFR_whole, TFR_ext, TFR_dir_wod) %>% 
  filter(Year >= 1751) %>%
  ggplot(aes(x = Year, y = TFR)) +
  geom_line(aes(colour = Dataset, linetype = Dataset), linewidth = 1.3)+ 
  scale_color_manual(values = c("#5B4792", "#621528", "#1A7761"))+
  scale_linetype_manual(values = c("22", "11", "solid")) +
  theme_graphs() 
# labs(title = "Total Fertility Rate in Sweden (1751-2022), 
#retrieved from a SOCSIM simulation and subsets of "extended" and direct" family trees without duplicates") 
ggsave(file="Graphs/Socsim_Exp2_TFR.jpeg", width=17, height=9, dpi=200)

# Life Expectancy at birth ----
# Calculate life expectancy at birth from asmr 1x1

# Estimate age-specific mortality rates 1x1 for the subset of direct ancestors and lateral kin without duplicates
asmr_ext_1 <- map_dfr(anc_kin_ext, ~ estimate_mortality_rates(opop = .x,
                                                              final_sim_year = 2022, #[Jan-Dec]
                                                              year_min = 1750, # Closed
                                                              year_max = 2023, # Open )
                                                              year_group = 1,
                                                              age_max_mort = 110, # Open )
                                                              age_group = 1), # [,)
                      .id = "Sim_id") 
save(asmr_ext_1, file = "Measures/asmr_ext_1.RData")

## In the first years of the subset of family trees, at almost all ages 
# rates have values of 0, Infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
# Hence, the life table cannot be computed and must be filter after some year
# where rates are available for the whole age range. 

# Load asmr 1x1 for the subset of "extended" family trees. 
load("Measures/asmr_ext_1.RData")

# Filter after year without warnings in filter(between(rn, 1, max(which(mx > 0))))
asmr_ext_1_filt <- asmr_ext_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  filter(Year > 1808)

# Compute life table from asmr 1x1 for the subset of "extended" family trees without duplicates
lt_ext <- lt_socsim(asmr_ext_1_filt)
save(lt_ext, file = "Measures/lt_ext.RData")

# Load asmr 1x1 for the subset of "direct" family trees
load("Measures/asmr_dir_wod_1.RData")
# Filter after year without warnings in filter(between(rn, 1, max(which(mx > 0))))
asmr_dir_wod_1_filt <- asmr_dir_wod_1 %>% 
  mutate(Year = as.numeric(str_extract(year, "\\d+"))) %>% 
  filter(Year > 1826)

# Compute life table from asmr 1x1 for the subset of "direct" family trees without duplicates
lt_dir_wod <- lt_socsim(asmr_dir_wod_1_filt)
save(lt_dir_wod, file = "Measures/lt_dir_wod.RData")

# Wrangle data for plotting

# Load age-specific mortality rates 1x1 for the whole single simulation (seed "13486"),
load("Measures/asmr_whole_1.RData")
# Load life table from asmr 1x1 for the whole single simulation (seed "13486"), used in 3_Compare_Ancestors.R
load("Measures/lt_whole.RData")

# Load life table from asmr 1x1 for the subset of "extended" family trees without duplicates
load("Measures/lt_ext.RData")

# Load life table from asmr 1x1 for the subset of "direct" family trees without duplicates
load("Measures/lt_dir_wod.RData")

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
lt_dir_wod2 <- lt_dir_wod %>%
  mutate(Year = as.numeric(str_extract(year, "\\d+")),
         Dataset = "Trees_direct",
         Rate = "e0") %>% 
  select(Year, ex, Dataset, Rate, sex, Age)

bind_rows(lt_whole2, lt_ext2, lt_dir_wod2) %>% 
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
ggsave(file="Graphs/socsim_Exp2_e0.jpeg", width=17, height=9, dpi=200)


#----------------------------------------------------------------------------------------------------
## Final plot combining TFR and e0 ----

## TFR and e0 (for females) from whole SOCSIM simulation and subsets of "extended" and direct" family trees 

yrs_plot2 <- c(1750, 1800, 1850, 1900, 1950, 2000)

bind_rows(TFR_whole %>% 
            rename(Estimate = TFR), 
          TFR_ext %>% 
            rename(Estimate = TFR), 
          TFR_dir_wod %>%           
            rename(Estimate = TFR)) %>%  
  bind_rows(lt_whole2 %>% 
              rename(Estimate = ex) %>% 
              filter(Age == 0),
            lt_ext2 %>% 
              rename(Estimate = ex) %>% 
              filter(Age == 0),
            lt_dir_wod2 %>% 
              rename(Estimate = ex) %>% 
              filter(Age == 0)) %>% 
  filter(sex == "female" & Year >= 1870) %>% 
  mutate(Rate = ifelse(Rate == "TFR", "Total Fertility Rate", "Life expectancy at birth"), 
         Rate = factor(Rate, levels = c("Total Fertility Rate", "Life expectancy at birth")), 
         Dataset = case_when(Dataset == "Trees_extended" ~ "Experiment 2 with lateral kin",
                             Dataset == "Trees_direct" ~ "Experiment 2 without lateral kin",
                             TRUE ~ "Whole simulation")) %>%
  ggplot(aes(x = Year, y = Estimate, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Rate, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), 
             aes(shape = Dataset), size = 6, fill = "#814352")+
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#8C7EB2", "#814352", "#1A7761"))+
  scale_shape_manual(values = c(18, 25, 46)) +
  scale_x_continuous(breaks = yrs_plot2)+
  theme_graphs()
# labs(title = "Total Fertility Rate and Life Expectancy at Birth in Sweden (1751-2020), retrieved 
# from SOCSIM microsimulation and subsets of "extended" and direct" family trees ") + 
# Save the plot
ggsave(file="Graphs/Final_Socsim_Exp2_TFR_e0.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Sex Ratio at Birth and Infant Mortality Rate ----

# Convert SOCSIM months to calendar years. 
asYr <- function(month, last_month, final_sim_year) {
  return(final_sim_year - trunc((last_month - month)/12))
}

## Define years if not set in the Global Environment
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

# Sex Ratio at Birth by year for "extended" family trees subset without duplicates
SRB_ext <- kin_ext %>% 
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
SRB_dir_wod <- kin_dir_wod %>% 
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
bind_rows(SRB_whole, SRB_ext, SRB_dir_wod) %>% 
  filter(Year >= 1751 & !is.nan(SRB) & !is.infinite(SRB)) %>% # No births after 2003??
  ggplot(aes(x = Year, y = SRB, group = Dataset, color = Dataset, linetype = Dataset))+
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#5B4792", "#621528", "#1A7761"))+
  scale_linetype_manual(values = c("22", "11", "solid")) +
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
Births_ext <- kin_ext %>% 
  mutate(Year = asYr(dob, last_month, final_sim_year)) %>% 
  filter(Year %in% year_range) %>% 
  count(Year) %>%
  mutate(Year = factor(Year, levels = year_range)) %>% 
  complete(Year, fill = list(n = 0))  %>%
  mutate(Year = as.numeric(as.character(Year)),
         Dataset = "Trees_extended", 
         Event = "Births")

# Births by year from subset of "direct" family trees without duplicates
Births_dir_wod <- kin_dir_wod %>% 
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
Deaths_0_ext <- kin_ext %>% 
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
Deaths_0_dir_wod <- kin_dir_wod %>% 
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
IMR <- bind_rows(Births_whole, Births_ext, Births_dir_wod, Deaths_0_whole, Deaths_0_ext, Deaths_0_dir_wod) %>%
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
bind_rows(SRB_whole, SRB_ext, SRB_dir_wod) %>% 
  select(Year, Dataset, SRB) %>% 
  full_join(IMR %>% select(Year, Dataset, IMR), by = c("Year", "Dataset")) %>% 
  pivot_longer(SRB:IMR, names_to = "Measure", values_to = "Value") %>% 
  # There can be infinite (x_Births_M/0_Births_F) 
  # and NaN (x_Births_M/0_Births_F or 0_Deaths/0_Births) values
  filter(Year >= 1870 & !is.infinite(Value) & !is.nan(Value)) %>%
  mutate(Measure = ifelse(Measure == "SRB", "Sex Ratio at Birth", "Infant Mortality Rate"), 
         Measure = factor(Measure, levels = c("Sex Ratio at Birth", "Infant Mortality Rate")), 
         Dataset = case_when(Dataset == "Trees_extended" ~ "Experiment 2 with lateral kin",
                             Dataset == "Trees_direct" ~ "Experiment 2 without lateral kin",
                             TRUE ~ "Whole simulation")) %>%
  ggplot(aes(x = Year, y = Value, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Measure, scales = "free") + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), 
             aes(shape = Dataset), size = 6, fill = "#814352")+
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#8C7EB2", "#814352", "#1A7761"))+
  theme_graphs() +
  scale_shape_manual(values = c(18, 25, 46)) +
  scale_x_continuous(breaks = yrs_plot2)+
  theme_graphs() +
  theme(axis.title.y = element_blank())
ggsave(file="Graphs/Socsim_Exp2_SRB_IMR.jpeg", width=17, height=9, dpi=200)

#----------------------------------------------------------------------------------------------------
## Births and Deaths counts by year -----

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
Deaths_ext <- kin_ext %>% 
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
Deaths_dir_wod <- kin_dir_wod %>% 
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
bind_rows(Births_whole, Births_ext, Births_dir_wod, Deaths_whole, Deaths_ext, Deaths_dir_wod) %>% 
  filter(Year >= 1870) %>%
  mutate(Dataset = case_when(Dataset == "Trees_extended" ~ "Experiment 2 with lateral kin",
                             Dataset == "Trees_direct" ~ "Experiment 2 without lateral kin",
                             TRUE ~ "Whole simulation")) %>%
  ggplot(aes(x = Year, y = n, group = Dataset, color = Dataset))+
  facet_wrap(. ~ Event) + 
  geom_point(data = . %>% filter(Year %in% yrs_plot2), aes(shape = Dataset), size = 7)+
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#5B4792", "#621528", "#1A7761"))+
  scale_shape_manual(values = c(15,17,46)) +
  theme_graphs() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(file="Graphs/Socsim_Exp2_Births_Deaths.jpeg", width=17, height=9, dpi=200)