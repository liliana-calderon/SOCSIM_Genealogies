#----------------------------------------------------------------------------------------------------
# Function to get pids of relevant kin up to the 4th degree of consanguinity of a given ego(s) 
# from SOCSIM microsimulation outputs

# To run the functions, .opop file must be set in the GlobalEnv

# Created by Liliana Calderon on 13-04-2022
# Last modified by Liliana Calderon on 17-03-2023

#------------------------------------------------------------------------------------------------------

## Function to replace zeros and empty vectors for NA 

# SOCSIM uses 0 (i.e., nothing) as the pid for no one. 
# I added the condition length(x)==0 to the original zna() function 
# as it returned logical (empty integer vectors) 
# for those members whose relative has already been replaced by NA
# This will not happen when tracing direct ancestors, but better to have a unique function

ze_na <- function(x) {
  x[x==0] <- NA
  x[length(x)==0] <- NA
  return(x) 
}


## Function to get pids of relevant kin up to the 4th degree of consanguinity from a given ego ----
# It takes as input a vector of ego's pids.


get_kin <- function(egos) {
  
  # Select pid of ego
  ego <- egos
  
  ### Find pid for relevant kin of ego.
  ### Each generation backwards is added to the right of the most recent generation
  
  ## 0-IL. Spouses
  w <- omar %>% filter(hpid == ego) %>% pull(wpid) %>% ze_na() # Wifes w marid
  h <- omar %>% filter(wpid == ego) %>% pull(hpid) %>% ze_na() # Husbands w marid
  
  
  ## 1. Parents
  m <- opop %>% filter(pid == ego) %>% pull(mom) %>% ze_na()
  f <- opop %>% filter(pid == ego) %>% pull(pop) %>% ze_na()
  
  ## 1. Children of ego 
  c <- opop %>% 
    filter(mom == ego | pop == ego) %>% 
    pull(pid) %>% 
    ze_na()
  
  ## 1-IL. Step-parents
  fw <- omar %>%
    filter(hpid == f & wpid != m) %>% # Father's wifes which are not ego's mother
    pull(wpid) %>%
    ze_na()
  mh <- omar %>%
    filter(wpid == m & hpid != f) %>% # Mother's husbands which are not ego's father
    pull(hpid) %>%
    ze_na()
  
  ## 1-IL. Step-children (children of all ego's spouses which are not from ego)
  sc <- opop %>%
    filter((mom %in% w & pop != ego) | (pop %in% h & mom != ego)) %>%
    pull(pid) %>%
    ze_na()
  
  ## 2. Grandparents
  mm <- opop %>% filter(pid == m) %>% pull(mom) %>% ze_na()
  mf <- opop %>% filter(pid == m) %>% pull(pop) %>% ze_na()
  fm <- opop %>% filter(pid == f) %>% pull(mom) %>% ze_na()
  ff <- opop %>% filter(pid == f) %>% pull(pop) %>% ze_na()
  
  ## 2-IL. Step-grand-parents. Parents of step-parents
  fwm <- opop %>% filter(pid %in% fw) %>% pull(mom) %>% ze_na()
  fwf <- opop %>% filter(pid %in% fw) %>% pull(pop) %>% ze_na()
  mhm <- opop %>% filter(pid %in% mh) %>% pull(mom) %>% ze_na()
  mhf <- opop %>% filter(pid %in% mh) %>% pull(pop) %>% ze_na()
    
  ## 2. Siblings (from both parents)
  z <- opop %>%
    filter(mom == m & pop == f & pid != ego) %>%
    pull(pid) %>%
    ze_na()
  
  ## 2. Half-Siblings (from only one parent)
  hz <- opop %>%
    filter((mom == m & pop != f) | (mom != m & pop == f)) %>%
    pull(pid) %>%
    ze_na()

  ## 2-IL. Step-siblings (children of step-parents)
  spc <- opop %>%
    filter((mom %in% fw & pop != f) | (pop %in% mh & mom != m)) %>%
    pull(pid) %>%
    ze_na()
  
  ## 2. Grand-children (Children of all ego's children)
  cc <- opop %>% 
    filter(mom %in% c | pop %in% c) %>% 
    pull(pid) %>% 
    ze_na()
  
  ## 2-IL. Step-grand-children (Children of ego's step-children)
  scc <- opop %>%
    filter((mom %in% sc) | (pop %in% sc)) %>%
    pull(pid) %>%
    ze_na()
  
  
  ## 3. Great-grandparents (1x).
  mmm <- opop %>% filter(pid == mm) %>% pull(mom) %>% ze_na()
  mmf <- opop %>% filter(pid == mm) %>% pull(pop) %>% ze_na()
  mfm <- opop %>% filter(pid == mf) %>% pull(mom) %>% ze_na()
  mff <- opop %>% filter(pid == mf) %>% pull(pop) %>% ze_na()
  fmm <- opop %>% filter(pid == fm) %>% pull(mom) %>% ze_na()
  fmf <- opop %>% filter(pid == fm) %>% pull(pop) %>% ze_na()
  ffm <- opop %>% filter(pid == ff) %>% pull(mom) %>% ze_na() 
  fff <- opop %>% filter(pid == ff) %>% pull(pop) %>% ze_na()
  
  ## 3. Aunts/uncles through grand-parents
  mz <- opop %>%   # From both maternal grandparents
    filter(mom == mm & pop == mf & pid != m) %>%
    pull(pid) %>%
    ze_na()
  mmc <- opop %>%  # From only one of the maternal grandparents
    filter(mom == mm & !pid %in% c(m, mz)) %>% 
    pull(pid) %>% 
    ze_na()
  mfc <- opop %>%  # From only one of the maternal grandparents
    filter(pop == mf & !pid %in% c(m, mz)) %>% 
    pull(pid) %>%
    ze_na()
  fz <- opop %>%   # From both Paternal grandparents
    filter(mom == fm & pop == ff & pid != f) %>%
    pull(pid) %>%
    ze_na()
  fmc <- opop %>%   # From only one of the paternal grandparents
    filter(mom == fm & !pid %in% c(f, fz)) %>% 
    pull(pid) %>% 
    ze_na()
  ffc <- opop %>%   # From only one of the paternal grandparents
    filter(pop == ff & !pid %in% c(f, fz)) %>% 
    pull(pid) %>% 
    ze_na()
  
  ## 3. Nieces/Nephews through sets of siblings
  zc <- opop %>% # From both parents siblings
    filter(mom %in% z | pop %in% z) %>%
    pull(pid) %>%
    ze_na()
  hzc <- opop %>% # From half-siblings
    filter(mom %in% hz | pop %in% hz) %>%
    pull(pid) %>%
    ze_na()

  ## 3. Great-grand-children. Children of grand-children
  ccc <- opop %>% 
    filter(mom %in% cc | pop %in% cc) %>% 
    pull(pid) %>% 
    ze_na()
  
  ## 4. Great-great-grandparents (2X)
  mmmm <- opop %>% filter(pid == mmm) %>% pull(mom) %>% ze_na()
  mmmf <- opop %>% filter(pid == mmm) %>% pull(pop) %>% ze_na()
  mmfm <- opop %>% filter(pid == mmf) %>% pull(mom) %>% ze_na()
  mmff <- opop %>% filter(pid == mmf) %>% pull(pop) %>% ze_na()
  mfmm <- opop %>% filter(pid == mfm) %>% pull(mom) %>% ze_na()
  mfmf <- opop %>% filter(pid == mfm) %>% pull(pop) %>% ze_na()
  mffm <- opop %>% filter(pid == mff) %>% pull(mom) %>% ze_na()
  mfff <- opop %>% filter(pid == mff) %>% pull(pop) %>% ze_na()
  
  fmmm <- opop %>% filter(pid == fmm) %>% pull(mom) %>% ze_na()
  fmmf <- opop %>% filter(pid == fmm) %>% pull(pop) %>% ze_na()
  fmfm <- opop %>% filter(pid == fmf) %>% pull(mom) %>% ze_na()
  fmff <- opop %>% filter(pid == fmf) %>% pull(pop) %>% ze_na()
  ffmm <- opop %>% filter(pid == ffm) %>% pull(mom) %>% ze_na()
  ffmf <- opop %>% filter(pid == ffm) %>% pull(pop) %>% ze_na()
  fffm <- opop %>% filter(pid == fff) %>% pull(mom) %>% ze_na()
  ffff <- opop %>% filter(pid == fff) %>% pull(pop) %>% ze_na()
  
  ## 4. Great-aunts/Great-uncles through great-grandparents (1X)
  
  mmz <- opop %>%   # From both maternal great-grandparents through maternal lineage
    filter(mom == mmm & pop == mmf & pid != mm) %>%
    pull(pid) %>%
    ze_na()
  
  mmhz <- opop %>%  # From only one maternal great-grandparent through maternal lineage, without distinction
    filter(mom == mmm | pop == mmf) %>% 
    filter(!pid %in% c(mm, mmz)) %>% 
    pull(pid) %>% 
    ze_na()
  
  mfz <- opop %>%   # From both paternal great-grandparents through maternal lineage
    filter(mom == mfm & pop == mff & pid != mf) %>%
    pull(pid) %>%
    ze_na()
  
  mfhz <- opop %>%  # From only one paternal great-grandparent through maternal lineage, without distinction
    filter(mom == mfm | pop == mff) %>% 
    filter(!pid %in% c(mf, mfz)) %>% 
    pull(pid) %>%
    ze_na()
  
  fmz <- opop %>%   # From both maternal great-grandparents through paternal lineage
    filter(mom == fmm & pop == fmf & pid != fm) %>%
    pull(pid) %>%
    ze_na()
  
  fmhz <- opop %>%   # From only one maternal great-grandparents through paternal lineage
    filter(mom == fmm | pop == fmf) %>% 
    filter(!pid %in% c(fm, fmz)) %>% 
    pull(pid) %>% 
    ze_na()
  
  ffz <- opop %>%   # From both paternal great-grandparents through paternal lineage
    filter(mom == ffm & pop == fff & pid != ff) %>%
    pull(pid) %>%
    ze_na()
  
  ffhz <- opop %>%   # From only one paternal great-grandparents through paternal lineage
    filter(mom == ffm | pop == fff) %>% 
    filter(!pid %in% c(ff, ffz)) %>% 
    pull(pid) %>% 
    ze_na()
  
  
  ## 4. Cousins from sets of aunts/uncles (defined by grandparents)
  mzc <- opop %>% # From both maternal grandparents
    filter(mom %in% mz | pop %in% mz) %>%
    pull(pid) %>%
    ze_na()
  mmcc <- opop %>% # From only one of the maternal grandparents
    filter(mom %in% mmc | pop %in% mmc) %>%
    pull(pid) %>%
    ze_na()
  mfcc <- opop %>% # From only one of the maternal grandparents
    filter(mom %in% mfc | pop %in% mfc) %>%
    pull(pid) %>%
    ze_na()
  fzc <- opop %>% # From both paternal grandparents
    filter(mom %in% fz | pop %in% fz) %>%
    pull(pid) %>%
    ze_na()
  fmcc <- opop %>% # From only one of the paternal grandparents
    filter(mom %in% fmc | pop %in% fmc) %>%
    pull(pid) %>%
    ze_na()
  ffcc <- opop %>% # From only one of the paternal grandparents
    filter(mom %in% ffc | pop %in% ffc) %>%
    pull(pid) %>%
    ze_na()
  
  ## Bind all ego's kin in a tibble and include their information from .opop
  
  kin_opop <- tibble(kin_type = c("ego",
                                  rep("w", length(w)),
                                  rep("h", length(h)),
                                  "m", "f",
                                  rep("c", length(c)),
                                  rep("fw", length(fw)),
                                  rep("mh", length(mh)),
                                  rep("sc", length(sc)),
                                  "mm", "mf", "fm", "ff",
                                  rep("fwm", length(fwm)),
                                  rep("fwf", length(fwf)),
                                  rep("mhm", length(mhm)),
                                  rep("mhf", length(mhf)),
                                  rep("z", length(z)), 
                                  rep("hz", length(hz)), 
                                  rep("spc", length(spc)),
                                  rep("cc", length(cc)),
                                  rep("scc", length(scc)),
                                  "mmm", "mmf", "mfm", "mff",
                                  "fmm", "fmf", "ffm", "fff",
                                  rep("mz", length(mz)), 
                                  rep("mmc", length(mmc)),
                                  rep("mfc", length(mfc)),
                                  rep("fz", length(fz)),
                                  rep("fmc", length(fmc)),
                                  rep("ffc", length(ffc)), 
                                  rep("zc", length(zc)), 
                                  rep("hzc", length(hzc)), 
                                  rep("ccc", length(ccc)), 
                                  "mmmm", "mmmf", "mmfm", "mmff", "mfmm", "mfmf", "mffm", "mfff", 
                                  "fmmm", "fmmf", "fmfm", "fmff", "ffmm", "ffmf", "fffm", "ffff", 
                                  rep("mmz", length(mmz)),
                                  rep("mmhz", length(mmhz)), 
                                  rep("mfz", length(mfz)), 
                                  rep("mfhz", length(mfhz)),
                                  rep("fmz", length(fmz)), 
                                  rep("fmhz", length(fmhz)), 
                                  rep("ffz", length(ffz)), 
                                  rep("ffhz", length(ffhz)), 
                                  rep("mzc", length(mzc)),
                                  rep("mmcc", length(mmcc)),
                                  rep("mfcc", length(mfcc)),
                                  rep("fzc", length(fzc)),
                                  rep("fmcc", length(fmcc)),
                                  rep("ffcc", length(ffcc))
                                  ), 
                     pid = c(ego,
                             w, h,
                             m, f, c,
                             fw, mh, sc,
                             mm, mf, fm, ff,
                             fwm, fwf, mhm, mhf,
                             z, hz,
                             spc, 
                             cc, 
                             scc,
                             mmm, mmf, mfm, mff,
                             fmm, fmf, ffm, fff,
                             mz, mmc, mfc, 
                             fz, fmc, ffc, 
                             zc, hzc, 
                             ccc,
                             mmmm, mmmf, mmfm, mmff, mfmm, mfmf, mffm, mfff, 
                             fmmm, fmmf, fmfm, fmff, ffmm, ffmf, fffm, ffff, 
                             mmz,mmhz, mfz, mfhz,
                             fmz, fmhz, ffz, ffhz, 
                             mzc, mmcc, mfcc, 
                             fzc, fmcc, ffcc),
                     ego_id = rep(ego, length(pid))) %>% 
    filter(!is.na(pid)) %>% 
    distinct(pid, .keep_all = TRUE)
  
  return(kin_opop)
  
}