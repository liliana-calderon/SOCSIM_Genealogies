#----------------------------------------------------------------------------------------------------
# Functions to Write SOCSIM fertility and mortality input rate files using HFD and HMD data
# retrieved using the HMDHFDplus package 
# and the Human Fertility COllection for the period not covered in HFD

## Write SOCSIM input rates from HMD and HFD using the HMDHFDplus package ----

# Created on 08-06-2022
# Last modified on 25-05-2023

#----------------------------------------------------------------------------------------------------
# Create a sub-folder called "rates" to save the rate files if it does not exist.
ifelse(!dir.exists("rates"), dir.create("rates"), FALSE)
# If the sub-folder name changes, 
# it must be also changed in the functions when opening the output file connection

#----------------------------------------------------------------------------------------------------
#### Write SOCSIM input fertility rates from HFD using the HMDHFDplus package ----

write_socsim_rates_HFD <- function(Country, HFD_username, HFD_password) {
  
  # Get ASFR from HFD
  asfr <- 
    readHFDweb(CNTRY = Country,
               item = "asfrRR", 
               username = HFD_username,
               password = HFD_password)
  
  # Wrangle data and compute monthly fertility rates
  ASFR <- 
    asfr %>% 
    select(-"OpenInterval") %>% 
    mutate(Age_up = Age + 1, # SOCSIM uses the upper age bound
           Month = 0, 
           ASFR_mo = ASFR/12) %>% 
    select(-ASFR)
  
  # Add rows with rates = 0 for ages 0-12 and 56-100
  ASFR <- 
    ASFR %>% 
    group_by(Year) %>% 
    group_split() %>% 
    map_df(~ add_row(.x,
                     Year = unique(.x$Year), 
                     Age = 0, Age_up = 12,  Month = 0, ASFR_mo = 0.0, 
                     .before = 1)) %>% 
    group_by(Year) %>% 
    group_split() %>% 
    map_df(~ add_row(.x, 
                     Year = unique(.x$Year), 
                     Age = 56, Age_up = 111, Month = 0, ASFR_mo = 0.0, 
                     .after = 45)) %>% 
    ungroup() %>% 
    select(-Age)
  
  # Extract the years available in HFD
  years <- ASFR %>% pull(Year) %>% unique()
  
  # Row numbers corresponding to sequence of years of age in ASFR
  rows_ageF <- ASFR %>% pull(Age_up) %>% unique() %>% seq_along()
  
  
  ## Write the fertility rate files for each year
  
  for(Year in years) {
    
    # Find the index of each year of the iteration
    n <- which(Year == years)
    n_row <- (n-1)*46 + rows_ageF
    
    # Open an output file connection
    outfilename <- file(paste0("rates/",Country,"fert",Year), "w") 
    # without ".txt" specification as the original files had no format. 
    
    # Include country and year of the corresponding rates
    cat(c("** Period (Monthly) Age-Specific Fertility Rates for", Country, "in", Year, "\n"), 
        file = outfilename)
    cat(c("* Retrieved from the Human Fertility Database, www.humanfertility.org", "\n"), 
        file = outfilename)
    cat(c("* Max Planck Institute for Demographic Research (Germany) and", "\n"), 
        file = outfilename)
    cat(c("* Vienna Institute of Demography (Austria)", "\n"), 
        file = outfilename)
    cat(c("* Data downloaded on ", format(Sys.time(), format= "%d %b %Y %X %Z"), 
          "using the HMDHFDplus R-Package", "\n"), 
        file = outfilename)
    cat(c("** NB: The original HFD annual rates have been converted into monthly rates", "\n"), 
        file = outfilename)
    cat(c("** The open age intervals (12- and 55+) are limited to one year [12-13) and [55-56)", "\n"), 
        file = outfilename)
    cat("\n", file = outfilename)
    
    # Print birth rates (single females)
    cat("birth", "1", "F", "single", "0", "\n", file = outfilename)
    for(i in n_row) {
      cat(c(as.matrix(ASFR)[i,-1], "\n"), file = outfilename) }
    cat("\n", file = outfilename)
    
    # Print birth rates (married females)
    cat("birth", "1", "F", "married", "0", "\n", file = outfilename)
    for(i in n_row) {
      cat(c(as.matrix(ASFR)[i,-1],"\n"), file = outfilename) }
    
    close(outfilename)
    
  }
  
}

#----------------------------------------------------------------------------------------------------
#### Write SOCSIM input mortality rates from HMD using the HMDHFDplus package ----

write_socsim_rates_HMD <- function(Country, HMD_username, HMD_password) {
  
  # Get female life tables from HMD
  ltf <- 
    readHMDweb(CNTRY = Country,
               item = "fltper_1x1",
               username = HMD_username,
               password = HMD_password)
  
  # Get male life tables from HMD
  ltm <- 
    readHMDweb(CNTRY = Country,
               item = "mltper_1x1",
               username = HMD_username,
               password = HMD_password)
  
  # Wrangle data and compute monthly mortality probabilities
  ASMP <- 
    ltf %>% 
    select(Year, Age, qx) %>% 
    left_join(ltm %>% select(Year, Age, qx), 
              by = c("Year","Age"), suffix = c("_F","_M")) %>% 
    mutate(qx_Fmo = ifelse(Age == 110, qx_F/12, 1-(1-qx_F)^(1/12)),
           qx_Mmo = ifelse(Age == 110, qx_M/12, 1-(1-qx_M)^(1/12)), 
           Age_up = Age + 1, # SOCSIM uses the upper age bound
           Month = 0) %>% 
    select(c(Year, Age_up, Month, qx_Fmo, qx_Mmo))
  
  # Extract the years available in HMD
  years <- ASMP %>% pull(Year) %>% unique()
  
  # Row numbers corresponding to sequence of years of age in ASMP
  rows_ageM <- ASMP %>% pull(Age_up) %>% unique() %>% seq_along()
  
  
  ## Write the mortality rate files for each year
  
  for(Year in years) {
    
    # Find the index of each year of the iteration
    n <- which(Year == years)
    n_row <- (n-1)*111 + rows_ageM
    
    # Open an output file connection
    outfilename <- file(paste0("rates/",Country,"mort",Year), "w") 
    # without ".txt" specification as the original files had no format. 
    
    # Include country and year of the corresponding probabilities
    cat(c("** Period (Monthly) Age-Specific Probabilities of Death for", Country, "in", Year, "\n"), 
        file = outfilename)
    cat(c("* Retrieved from the Human Mortality Database (Life tables), www.mortality.org.", "\n"), 
        file = outfilename)
    cat(c("* Max Planck Institute for Demographic Research (Germany),", "\n"), 
        file = outfilename)
    cat(c("* University of California, Berkeley (USA),", "\n"), 
        file = outfilename)
    cat(c("* and French Institute for Demographic Studies (France)", "\n"), 
        file = outfilename)
    cat(c("* Data downloaded on ", format(Sys.time(), format= "%d %b %Y %X %Z"), 
          "using the HMDHFDplus R-Package", "\n"), 
        file = outfilename)
    cat(c("** NB: The original HMD annual probabilities have been converted into monthly probabilities", "\n"), 
        file = outfilename)
    cat(c("** The open age interval (110+) is limited to one year [110-111)", "\n"), 
        file = outfilename)
    cat("\n", file = outfilename)
    
    # Print mortality probabilities (single females)
    cat("death", "1", "F", "single", "\n", file = outfilename)
    for(i in n_row) {
      cat(c(as.matrix(ASMP)[i,-c(1,5)], "\n"), file = outfilename) }
    cat("\n", file = outfilename)
    
    # Print mortality probabilities (single males)
    cat("death", "1", "M", "single", "\n", file = outfilename)
    for(i in n_row) {
      cat(c(as.matrix(ASMP)[i,-c(1,4)], "\n"), file = outfilename) }
    
    close(outfilename)
    
  }
}

#----------------------------------------------------------------------------------------------------
#### Write SOCSIM input fertility rates from HFC for period not covered by HFD ----

write_socsim_rates_HFC <- function(Country, RefCode, Year_Max) {
  
  ## Download ASFR from HFC
  HFC <- read_csv(paste0("https://www.fertilitydata.org/File/GetFile/Country/",Country,"/",Country,"_ASFRstand_TOT.txt"),
                  col_names = T, show_col_types = F)
  
  # Wrangling HFC data
  HFC <- HFC %>% 
    filter(RefCode %in% "SWE_02" & Year2 <= Year_Max) %>% 
    select(Year1, Age, ASFR) %>% # Select useful columns
    mutate(ASFR = ifelse(ASFR == '.', "0", ASFR), 
           ASFR = as.numeric(ASFR),
           Age_up = Age + 1, # SOCSIM uses the upper age bound
           Month = 0, 
           ASFR_mo = ASFR/12) %>%
    select(-ASFR)
  
  # Add rows with rates = 0 for ages [0-12] and [51-111] to keep the same format than HFD
  
  ASFR <- 
    HFC %>% 
    group_by(Year1) %>% 
    group_split() %>% 
    map_df(~ add_row(.x,
                     Year1 = unique(.x$Year1), 
                     Age = 0, Age_up = 12,  Month = 0, ASFR_mo = 0.0, 
                     .before = 1)) %>% 
    group_by(Year1) %>% 
    group_split() %>% 
    map_df(~ add_row(.x,
                     Year1 = unique(.x$Year1), 
                     Age = 12, Age_up = 13,  Month = 0, ASFR_mo = 0.0, 
                     .after = 1)) %>% 
    group_by(Year1) %>% 
    group_split() %>% 
    map_df(~ add_row(.x,
                     Year1 = unique(.x$Year1), 
                     Age = 13, Age_up = 14,  Month = 0, ASFR_mo = 0.0, 
                     .after = 2)) %>% 
    group_by(Year1) %>% 
    group_split() %>% 
    map_df(~ add_row(.x, 
                     Year1 = unique(.x$Year1), 
                     Age = 51, Age_up = 52, Month = 0, ASFR_mo = 0.0,
                     .after = 40)) %>%
    group_by(Year1) %>% 
    group_split() %>% 
    map_df(~ add_row(.x, 
                     Year1 = unique(.x$Year1), 
                     Age = 52, Age_up = 53, Month = 0, ASFR_mo = 0.0, 
                     .after = 41)) %>%
    group_by(Year1) %>% 
    group_split() %>% 
    map_df(~ add_row(.x, 
                     Year1 = unique(.x$Year1), 
                     Age = 53, Age_up = 54, Month = 0, ASFR_mo = 0.0,
                     .after = 42)) %>%
    group_by(Year1) %>% 
    group_split() %>% 
    map_df(~ add_row(.x, 
                     Year1 = unique(.x$Year1), 
                     Age = 54, Age_up = 55, Month = 0, ASFR_mo = 0.0, 
                     .after = 43)) %>%
    group_by(Year1) %>% 
    group_split() %>% 
    map_df(~ add_row(.x, 
                     Year1 = unique(.x$Year1), 
                     Age = 55, Age_up = 56, Month = 0, ASFR_mo = 0.0, 
                     .after = 44)) %>%
    group_by(Year1) %>% 
    group_split() %>% 
    map_df(~ add_row(.x, 
                     Year1 = unique(.x$Year1), 
                     Age = 56, Age_up = 111, Month = 0, ASFR_mo = 0.0, 
                     .after = 45)) %>%
    ungroup() %>%
    # Add years that hold for each set of rates
    mutate(Year2 = Year1 + 1, 
           Year3 = Year1 + 2, 
           Year4 = Year1 + 3,
           Year5 = Year1 + 4) %>% 
    select(Year1, Year2, Year3, Year4, Year5, Age_up, Month, ASFR_mo) %>%
    # Repeat rates for year groups for each calendar year
    pivot_longer(cols = c(Year1:Year5), names_to = "Delete", values_to = "Year") %>% 
    select(Year, Age_up, Month, ASFR_mo) %>% 
    arrange(Year, Age_up)
  
  # Extract the years available in HFC
  years <- ASFR %>% pull(Year) %>% unique()
  
  # Row numbers corresponding to sequence of years of age in ASFR
  rows_ageF <- ASFR %>% pull(Age_up) %>% unique() %>% seq_along()
  
  ## Write the fertility rate files for each year
  
  for(Year in years) {
    
    # Find the index of each year of the iteration
    n <- which(Year == years)
    n_row <- (n-1)*46 + rows_ageF
    
    # Open an output file connection
    outfilename <- file(paste0("rates/",Country,"fert",Year), "w") 
    # without ".txt" specification as the original files had no format. 
    
    # Include country and year of the corresponding rates
    cat(c("** Period (Monthly) Age-Specific Fertility Rates for", Country, "in", Year, "\n"), 
        file = outfilename)
    cat(c("* Retrieved from the Human Fertility Collection, www.fertilitydata.org", "\n"), 
        file = outfilename)
    cat(c("* Max Planck Institute for Demographic Research (Germany) and", "\n"), 
        file = outfilename)
    cat(c("* Vienna Institute of Demography (Austria)", "\n"), 
        file = outfilename)
    cat(c("* Data downloaded on ", format(Sys.time(), format= "%d %b %Y %X %Z"), "\n"), 
        file = outfilename)
    cat(c("** NB: The original HFC annual rates have been converted into monthly rates", "\n"), 
        file = outfilename)
    cat(c("** The open age intervals (14- and 49+) are limited to one year [14-15) and [49-50)", "\n"), 
        file = outfilename)
    cat(c("** The ages below and above are set to 0", "\n"), file = outfilename)
    cat("\n", file = outfilename)
    
    # Print birth rates (single females)
    cat("birth", "1", "F", "single", "0", "\n", file = outfilename)
    for(i in n_row) {
      cat(c(as.matrix(ASFR)[i,-1], "\n"), file = outfilename) }
    cat("\n", file = outfilename)
    
    # Print birth rates (married females)
    cat("birth", "1", "F", "married", "0", "\n", file = outfilename)
    for(i in n_row) {
      cat(c(as.matrix(ASFR)[i,-1],"\n"), file = outfilename) }
    
    close(outfilename)
    
  }
  
}