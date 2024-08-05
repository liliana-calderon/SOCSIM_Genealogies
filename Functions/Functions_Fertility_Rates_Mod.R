#----------------------------------------------------------------------------------------------------
# Functions to estimate age-specific fertility from SOCSIM output ----

# Created on 18-01-2022
# Last modified on 22-04-2024

#----------------------------------------------------------------------------------------------------
## Notes ----

# We use mid-year population size by sex and age as an estimate of person-years lived during the year.
# Hence, the population includes all those who were born before the 1st of July of a given year 
# and die in or after July of that year or are still alive at the end of the simulation.
# This calculation could still be corrected to obtain a more precise estimate of exposures-to-risk
# e.g., HFD/HMD exposures-to-risk are based on annual population estimates (1st January of years t and t+1)
# with small corrections to account for the timing of deaths during the interval (Upper and Lower Lexis triangles)
# as well as the variations in the distribution of births by months.

#----------------------------------------------------------------------------------------------------
## Functions to convert SOCSIM months into calendar time ----

# Convert SOCSIM months to calendar years. 
asYr <- function(month, last_month, final_sim_year) {
  return(final_sim_year - trunc((last_month - month)/12))
}

## Returns simulation month corresponding to July of a given real calendar year
jul <- function(year, last_month, final_sim_year){
  return((last_month-5) - (final_sim_year - year)*12)
}

#----------------------------------------------------------------------------------------------------
#### Estimate SOCSIM period age-specific fertility rates ----
# This is a modified version of the function included in the rsocsim package
# It allows to handle the errors emerging from intentional duplicates in the data
# specially when left_joining the opop data frame to add the mothers' birth year
# which are duplicated in the opop file

estimate_fertility_rates_mod <- function(opop, final_sim_year, year_min, year_max, year_group, age_min_fert, age_max_fert, age_group) {
  
  last_month <- max(opop$dod) # Change to dod as max(dob) in most subsets and simulations is not the last simulated month
  
  # Year range and breaks
  year_range <- year_min:(year_max-1)
  year_breaks <- seq(year_min, year_max, by = year_group)
  
  # Age groups of fertility rates
  age_breaks_fert <- seq(age_min_fert, age_max_fert, by = age_group)
  
  # Calculate year of birth
  opop2 <- opop %>%
    mutate(birth_year = asYr(dob, last_month, final_sim_year))
  
  # 1. Numerator - births by women of given age group
  numerator <- yearly_birth_by_age_socsim(opop = opop2, 
                                          year_range = year_range, 
                                          age_breaks_fert = age_breaks_fert)
  
  # 2. Denominator - women in reproductive years (1st July)
  denom <- lapply(year_range, get_women_reproductive_age_socsim,
                  opop = opop2,
                  final_sim_year = final_sim_year,
                  age_breaks_fert = age_breaks_fert)
  
  denominator <- data.frame(do.call(rbind, denom))
  
  # 3. Rate
  asfr <- bind_cols(denominator %>% rename(deno = n),
                    numerator %>% select(nume = n)) %>%
    mutate(socsim = nume/deno,
           year_gr = cut(year, breaks = year_breaks, 
                         include.lowest = F, right = F, ordered_results = T)) %>%
    group_by(year = year_gr, age = agegr) %>%
    summarise(socsim = mean(socsim)) %>%
    ungroup
  
  return(asfr)
  
}

#----------------------------------------------------------------------------------------------------
#### Yearly Births by Age SOCSIM ----

yearly_birth_by_age_socsim <- function(opop, year_range, age_breaks_fert) {
  
  last_month <- max(opop$dod) # Change to dod as max(dob) in most subsets and simulations is not the last simulated month
  
  out <- opop %>% 
    left_join(opop %>% 
                # Add distinct function to avoid joining duplicated values (multiple = "first no longer works)
                distinct(pid, .keep_all = T) %>% 
                select(mom = pid, mother_birth = birth_year),
              by = "mom") %>% 
    select(birth_year, mother_birth) %>% 
    filter(birth_year %in% year_range) %>% 
    mutate(birth_year_factor = factor(birth_year, levels = year_range),
           mother_age = birth_year - mother_birth,
           mother_agegr_factor = cut(mother_age, breaks = age_breaks_fert, 
                                     include.lowest = F, right = F, ordered_results = T)) %>%
    filter(!is.na(mother_agegr_factor)) %>% 
    count(birth_year = birth_year_factor, mother_agegr = mother_agegr_factor) %>% 
    complete(birth_year, mother_agegr, fill = list(n = 0)) %>% 
    select(year = birth_year, agegr = mother_agegr, n) %>% 
    arrange(year, agegr)
  
  return(out)
  
}

#----------------------------------------------------------------------------------------------------
#### Get Women Reproductive Age SOCSIM (Female mid-year population) ----
# Return women by age alive on the 1st of July of a given year, including right-censored

get_women_reproductive_age_socsim <- function(opop, final_sim_year, year, age_breaks_fert) {
  
  last_month <- max(opop$dod) # Change to dod as max(dob) in most subsets and simulations is not the last simulated month
  opop$census <- year
  
  out <- opop %>% 
    mutate(dod2 = ifelse(dod == 0, 999999999, dod)) %>%
    filter(fem == 1 & dob < jul(year, last_month, final_sim_year) & dod2 >= jul(year, last_month, final_sim_year)) %>%
    mutate(age_at_census = trunc((jul(census, last_month, final_sim_year)-dob)/12),
           agegr_at_census = cut(age_at_census, breaks = age_breaks_fert, 
                                 include.lowest = F, right = F, ordered_results = T)) %>% 
    filter(!is.na(agegr_at_census)) %>% 
    count(agegr_at_census, census) %>% 
    complete(agegr_at_census, census, fill = list(n = 0)) %>% 
    select(year = census, agegr = agegr_at_census, n) %>% 
    arrange(year, agegr)
  
  return(out)
}