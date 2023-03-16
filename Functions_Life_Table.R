#----------------------------------------------------------------------------------------------------
# Function to calculate life table from SOCSIM output  ----

# Created by Liliana Calderon on 15-03-2023
# NB: Currently, this function only works with the asmr df calculated with the get_asmr_socsim()
# This code was inspired by Tim Riffe's BSSD2021Module2 code to calculate life tables
# https://github.com/timriffe/BSSD2021Module2/blob/master/02_tuesday/02_tuesday.Rmd 
# Last modified by Liliana Calderon on 16-03-2023
#----------------------------------------------------------------------------------------------------

lt_socsim <- function(asmr_socsim){
  
  lt_socsim <- asmr_socsim %>% 
    mutate(Age = as.numeric(str_extract(age, "\\d+"))) %>% 
    rename(mx = socsim) %>% 
    group_by(year, sex) %>% 
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
  
  return(lt_socsim)
  
}
