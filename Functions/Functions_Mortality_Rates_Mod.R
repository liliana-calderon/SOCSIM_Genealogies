#----------------------------------------------------------------------------------------------------
# Functions to estimate age-specific mortality from SOCSIM output ----

# Created on 22-04-2024
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

# Due to the limited population size in a microsimulation (especially, of survivors at old ages)
# sometimes few or no individuals of a specific age are alive at a exact point in time (here, 1st July). 
# Hence, it is possible to obtain rates higher than 1, equal to 0 (0_Events/Pop), 
# infinite (Events/0_Pop) and NaN (0_Events/0_Pop) values.

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
#### Get SOCSIM period ASMR ----

estimate_mortality_rates_mod <- function(opop, final_sim_year, year_min, year_max, year_group, age_max_mort, age_group) {
  
  last_month <- max(opop$dod) # Change to dod as minimum dob in this experiment is 18 years ago (age of genealogists)
  
  # Year range and breaks 
  year_range <- year_min:(year_max-1)
  year_breaks <- seq(year_min, year_max, by = year_group)
  
  # Age groups of mortality rates
  if(age_group == 1){
    age_breaks_mort <- c(0, seq(age_group, age_max_mort, by = age_group))
  } else{ age_breaks_mort <- c(0, 1, seq(age_group, age_max_mort, by = age_group))
  }
  
  # Age levels to complete the census (denominator)
  age_levels_census <- seq(0, age_max_mort-1, by = 1)
  
  # Calculate age at death and year of death
  opop2 <- opop %>% 
    rename(sex = fem) %>% 
    mutate(age_death_months = ifelse(dod == 0,NA,dod-dob), 
           age_death = trunc(age_death_months/12), 
           age_death_g = cut(age_death, breaks = age_breaks_mort, 
                             include.lowest = F, right = F, ordered_results = T),
           death_year = ifelse(dod == 0, NA, asYr(dod, last_month, final_sim_year)))
  
  # 1. Numerator - death counts by sex and age
  numerator <- opop2 %>% 
    filter(dod != 0 & death_year %in% year_range & !is.na(age_death_g)) %>%
    count(death_year, sex, age_death_g) %>% 
    rename(year = death_year, age = age_death_g) %>% 
    arrange(year, sex, age) %>% 
    mutate(year = factor(year, levels = year_range),
           sex = factor(sex, levels = c("0","1"))) %>% 
    complete(year, sex, age, fill = list(n = 0))  %>%
    mutate(year = as.numeric(as.character(year)))
  
  # 2. Denominator - population by sex and age (1st July)
  opop2_subset <- opop2 %>% 
    select(pid, dob, dod, sex)
  
  yearly_pop_age_sex <- lapply(year_range, census_socsim, 
                               opop = opop2_subset, 
                               final_sim_year = final_sim_year,
                               age_levels_census = age_levels_census)  
  
  denominator <- data.frame(do.call(rbind, yearly_pop_age_sex)) %>% 
    rename(year = census) %>% 
    mutate(age_at_census = as.numeric(as.character(age_at_census)),
           age = cut(age_at_census, breaks = age_breaks_mort, 
                     include.lowest = F, right = F, ordered_results = T)) %>%
    filter(!is.na(age)) %>% 
    group_by(year, sex, age) %>% 
    summarise(n = sum(n)) %>% 
    arrange(year, sex, age) %>% 
    ungroup()
  
  # 3. Rate
  asmr <- numerator %>% 
    full_join(denominator, by = c("year", "sex", "age"), 
              suffix = c("_num", "_den")) %>% 
    mutate(socsim = n_num / n_den,
           year_gr = cut(year, breaks = year_breaks, 
                         include.lowest = F, right = F, ordered_results = T),
           sex = ifelse(sex == 0, "male", "female")) %>%
    group_by(year = year_gr, sex, age) %>%
    summarise(socsim = mean(socsim)) %>%
    ungroup()
  
  return(asmr)  
  
}

#----------------------------------------------------------------------------------------------------
#### Census SOCSIM (Mid-year population, 1st July) ----
# Returns population by sex and age alive on the 1st of July of a given year, including right-censored

census_socsim <- function(opop, year, final_sim_year, age_levels_census) {
  
  last_month <- max(opop$dod) # Change to dod as minimum dob in this experiment is 18 years ago (age of genealogists)
  opop$census <- year 
  
  out <- opop %>% 
    mutate(dod2 = ifelse(dod == 0, 999999999, dod)) %>%
    filter(dob < jul(year, last_month, final_sim_year) & dod2 >= jul(year, last_month, final_sim_year)) %>% 
    mutate(age_at_census = trunc((jul(census, last_month, final_sim_year)-dob)/12)) %>%
    count(sex, age_at_census, census) %>% 
    mutate(sex = factor(sex, levels = c("0","1")),
           age_at_census = factor(age_at_census, levels = age_levels_census)) %>%
    complete(sex, age_at_census, census, fill = list(n = 0))
  
  return(out)
}