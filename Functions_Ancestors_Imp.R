#----------------------------------------------------------------------------------------------------
# Function to get pids of direct ancestors up to the 9th generation of a given ego(s)
# from SOCSIM microsimulation outputs

# To run the functions, .opop file must be set in the GlobalEnv

# Created by Liliana Calderon on 23-09-2022
# Last modified by Liliana Calderon on 31-05-2023
#----------------------------------------------------------------------------------------------------

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


## Function to get pids of direct ancestors from a given ego ----
# It takes as input a vector of ego's pids

get_ancestors <- function(egos) {
  
  # 1. Select pid of ego
  ego <- egos
  
  ### Find pid for relevant kin of ego.
  # NB: Each generation backwards is added to the left of the most recent generation

  ## 2. Parents (Mothers and fathers). 2
  m <- opop %>% filter(pid == ego) %>% pull(mom) %>% ze_na()
  f <- opop %>% filter(pid == ego) %>% pull(pop) %>% ze_na()
  
  ## 3. Grandparents. 4
  gm <- opop %>% filter(pid %in% c(m, f)) %>% pull(mom) %>% ze_na()
  gf <- opop %>% filter(pid %in% c(m, f)) %>% pull(pop) %>% ze_na()
  
  ## 4. Great-grandparents (1x). 8
  ggm <- opop %>% filter(pid %in% c(gm, gf)) %>% pull(mom) %>% ze_na()
  ggf <- opop %>% filter(pid %in% c(gm, gf)) %>% pull(pop) %>% ze_na()
  
  ## 5. Great-great-grandparents (2x). 16
  gggm <- opop %>% filter(pid %in% c(ggm, ggf)) %>% pull(mom) %>% ze_na()
  gggf <- opop %>% filter(pid %in% c(ggm, ggf)) %>% pull(pop) %>% ze_na()
  
  ## 6. Great-great-great-grandparents (3x). 32
  ggggm <- opop %>% filter(pid %in% c(gggm, gggf)) %>% pull(mom) %>% ze_na()
  ggggf <- opop %>% filter(pid %in% c(gggm, gggf)) %>% pull(pop) %>% ze_na()
  
  ## 7. Great-great-great-great-grandparents (4x). 64
  gggggm <- opop %>% filter(pid %in% c(ggggm, ggggf)) %>% pull(mom) %>% ze_na()
  gggggf <- opop %>% filter(pid %in% c(ggggm, ggggf)) %>% pull(pop) %>% ze_na()

  ## 8. Great-great-great-great-great-grandparents (5x). 128
  ggggggm <- opop %>% filter(pid %in% c(gggggm, gggggf)) %>% pull(mom) %>% ze_na()
  ggggggf <- opop %>% filter(pid %in% c(gggggm, gggggf)) %>% pull(pop) %>% ze_na()

  ## 9. Great-great-great-great-great-great-grandparents (6x). 256
  gggggggm <- opop %>% filter(pid %in% c(ggggggm, ggggggf)) %>% pull(mom) %>% ze_na()
  gggggggf <- opop %>% filter(pid %in% c(ggggggm, ggggggf)) %>% pull(pop) %>% ze_na()
  
 
  ## Bind all ego's direct ancestors in a tibble and include their information from .opop
  
  ancestors_opop <- tibble(kin_type = c("ego", # 1
                                        "mother", # 2
                                        "father", # 2
                                        rep("gmother", length(gm)), # 3 
                                        rep("gfather", length(gf)), # 3
                                        rep("ggmother", length(ggm)), # 4
                                        rep("ggfather", length(ggf)), # 4
                                        rep("gggmother", length(gggm)), # 5
                                        rep("gggfather", length(gggf)), # 5
                                        rep("ggggmother", length(ggggm)), # 6
                                        rep("ggggfather", length(ggggf)),# 6
                                        rep("gggggmother", length(gggggm)), # 7
                                        rep("gggggfather", length(gggggf)), # 7
                                        rep("ggggggmother", length(ggggggm)), # 8
                                        rep("ggggggfather", length(ggggggf)), # 8
                                        rep("gggggggmother", length(gggggggm)), # 9
                                        rep("gggggggfather", length(gggggggf))), # 9
                           pid = c(ego, # 1 
                                   m, f, # 2 
                                   gm, gf, # 3 
                                   ggm, ggf, # 4
                                   gggm, gggf, # 5
                                   ggggm, ggggf, # 6
                                   gggggm, gggggf, # 7
                                   ggggggm, ggggggf, # 8
                                   gggggggm, gggggggf # 9
                                   ),
  ego_id = rep(ego, length(pid))) %>% 
    filter(!is.na(pid)) %>% 
    distinct(pid, .keep_all = TRUE)

  return(ancestors_opop)
  
}