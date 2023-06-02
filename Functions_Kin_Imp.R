#----------------------------------------------------------------------------------------------------
# Function to get pids of relevant kin up to the 4th degree of consanguinity of a given ego(s) 
# from SOCSIM microsimulation outputs

# To run the functions, .opop file must be set in the GlobalEnv

# Created by Liliana Calderon on 13-04-2022
# Last modified by Liliana Calderon on 02-06-2023
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
  ## Number refers to the degree of consanguinity or affinity. Ego here is 0
  
  ## 0-IL. Spouses
  
  # w <- omar %>% filter(hpid == ego) %>% pull(wpid) %>% ze_na() # Wifes with marid
  # h <- omar %>% filter(wpid == ego) %>% pull(hpid) %>% ze_na() # Husbands w marid
  
  ## 1. Parents (Mothers and fathers).
  m <- opop %>% filter(pid == ego) %>% pull(mom) %>% ze_na()
  f <- opop %>% filter(pid == ego) %>% pull(pop) %>% ze_na()
  
  ## 1. Children of ego 
  c <- opop %>% filter(mom == ego | pop == ego) %>% pull(pid) %>% ze_na()

  ## 2. Grandparents
  gm <- opop %>% filter(pid %in% c(m, f)) %>% pull(mom) %>% ze_na()
  gf <- opop %>% filter(pid %in% c(m, f)) %>% pull(pop) %>% ze_na()

  ## 2. Siblings
  z <- opop %>% filter(mom == m & pop == f & pid != ego) %>% pull(pid) %>% ze_na()

  ## 2. Half siblings
  hz <- opop %>% filter((mom == m & pop != f) | (mom != m & pop == f)) %>%
    pull(pid) %>%  ze_na()
  
  ## 2. Grandchildren (Children of ego's children)
  gc <- opop %>% filter(mom %in% c | pop %in% c) %>% pull(pid) %>% ze_na()
  
  ## 3. Great-grandparents
  ggm <- opop %>% filter(pid %in% c(gm, gf)) %>% pull(mom) %>% ze_na()
  ggf <- opop %>% filter(pid %in% c(gm, gf)) %>% pull(pop) %>% ze_na()
  
  ## 3. Aunts/uncles (grand-parents' children)
  au <- opop %>% filter(mom %in% gm | pop %in% gf) %>% 
    filter(!pid %in% c(m, f)) %>% pull(pid) %>% ze_na()
  
  ## 3. Nieces/Nephews (siblings' children)
  ni <- opop %>% filter(mom %in% z | pop %in% z) %>% pull(pid) %>% ze_na()
  
  ## 3. Great-grand-children. (Grand-children's children)
  ggc <- opop %>% filter(mom %in% gc | pop %in% gc) %>% pull(pid) %>% ze_na()
  
  ## 4. Great-great-grandparents
  gggm <- opop %>% filter(pid %in% c(ggm, ggf)) %>% pull(mom) %>% ze_na()
  gggf <- opop %>% filter(pid %in% c(ggm, ggf)) %>% pull(pop) %>% ze_na()

  ## 4. Great-aunts/uncles (Great-grandparents' children)
  gau <- opop %>%   
    filter(mom %in% ggm | pop %in% ggf) %>%  
    filter(!pid %in% c(gm, gf)) %>% pull(pid) %>% ze_na()
  
  ## 4. Cousins (Aunts/uncles' children)
  k <- opop %>%   
    filter(mom %in% au | pop %in% au) %>% 
    filter(pid != ego) %>% pull(pid) %>% ze_na()
  
## In-law relatives
  
  ## 1-IL. Step-parents
  ## Father's wifes which are not ego's mother
  # sm <- omar %>% filter(hpid == f & wpid != m) %>% pull(wpid) %>% ze_na()
  ## Mother's husbands which are not ego's father
  # sf <- omar %>% filter(wpid == m & hpid != f) %>% pull(hpid) %>% ze_na()
  
  ## 1-IL. Step-children (children of all ego's spouses which are not from ego)
  # sc <- opop %>%
  #   filter((mom %in% w & pop != ego) | (pop %in% h & mom != ego)) %>%
  #   pull(pid) %>%
  #   ze_na()
  
  ## 2-IL. Step-grand-parents.
  # sgm <- opop %>% filter(pid %in% c(sm, sf)) %>% pull(mom) %>% ze_na()
  # sgf <- opop %>% filter(pid %in% c(sm, sf)) %>% pull(pop) %>% ze_na()
  
  ## 2-IL. Step-siblings (children of step-parents)
  # sz <- opop %>%
  #   filter((mom %in% sm & pop != f) | (pop %in% sf & mom != m)) %>%
  #   pull(pid) %>%
  #   ze_na()
  
  ## 2-IL. Step-grandchildren (Children of ego's step-children)
  # sgc <- opop %>% filter(mom %in% sc | pop %in% sc) %>% pull(pid) %>% 
  #   ze_na()
  
  ## 1-IL. Parents-in-law
  # ml <- opop %>% filter(pid %in% c(w, h)) %>% pull(mom) %>% ze_na()
  # fl <- opop %>% filter(pid %in% c(w, h)) %>% pull(pop) %>% ze_na()
  
  ## 1-IL. Children-in-law ??
  
  ## 2-IL. Siblings-in-law ??
  
  
  ## Bind all ego's kin in a tibble and include their information from .opop
  
  kin_opop <- tibble(kin_type = c("ego",
                                  # rep("wife", length(w)),
                                  # rep("husband", length(h)),
                                  "mother", 
                                  "father",
                                  rep("children", length(c)),
                                  rep("gmother", length(gm)),
                                  rep("gfather", length(gf)),
                                  rep("siblings", length(z)), 
                                  rep("half_siblings", length(hz)), 
                                  rep("gchildren", length(gc)), 
                                  rep("ggmother", length(ggm)),
                                  rep("ggfather", length(ggf)),
                                  rep("aunts_uncles", length(au)),
                                  rep("niblings", length(ni)),
                                  rep("ggchildren", length(ggc)),
                                  rep("gggmother", length(gggm)),
                                  rep("gggfather", length(gggf)), 
                                  rep("gaunts_guncles", length(gau)), 
                                  rep("cousins", length(k))
                                  ## step-relatives (up to 2nd degree of affinity)
                                  # rep("step_mother", length(sm)),
                                  # rep("step_father", length(sf)),
                                  # rep("step_children", length(sc)),
                                  # rep("step_gmother", length(sgm)),
                                  # rep("step_gfather", length(sgf)),
                                  # rep("step_siblings", length(sz)),
                                  # rep("step_gchildren", length(sgc)),
                                  # rep("mother_in_law", length(ml)),
                                  # rep("father_in_law", length(fl))
                                  ), 
                     pid = c(ego,
                             # w, h,
                             m, f, 
                             c,
                             gm, gf,
                             z, hz,
                             gc, 
                             ggm, ggf, 
                             au, ni, 
                             ggc, 
                             gggm, gggf, 
                             gau, 
                             k
                             # sm, sf, sc,
                             # sgm, sgf,
                             # sz, sgc, 
                             # ml, fl
                             ),
                     ego_id = rep(ego, length(pid))) %>% 
    filter(!is.na(pid)) %>% 
    distinct(pid, .keep_all = TRUE)
  
  return(kin_opop)
  
}