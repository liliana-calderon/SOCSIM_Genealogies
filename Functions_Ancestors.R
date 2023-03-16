#----------------------------------------------------------------------------------------------------
# Functions to trace direct ancestors of a given ego(s) 

# To run the functions, .opop file (and omar) must be set in the GlobalEnv

# Created by Liliana Calderon on 23-09-2022
# Last modified by Liliana Calderon on 14-03-2023
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

## Function to get pids and opop information of relevant kin from a given ego ----

get_ancestors <- function(egos) {
  
  # 1. Select pid of ego
  ego <- egos
  
  ### Find pid for relevant kin of ego
  
  ## 2. Parents
  m <- opop %>% filter(pid == ego) %>% pull(mom) %>% ze_na()
  f <- opop %>% filter(pid == ego) %>% pull(pop) %>% ze_na()
  
  # NB: Each generation backwards is added to the right of the most recent generation
  
  ## 3. Maternal grandparents
  mm <- opop %>% filter(pid == m) %>% pull(mom) %>% ze_na()
  mf <- opop %>% filter(pid == m) %>% pull(pop) %>% ze_na()
  
  ## 3. Paternal grandparents
  fm <- opop %>% filter(pid == f) %>% pull(mom) %>% ze_na()
  ff <- opop %>% filter(pid == f) %>% pull(pop) %>% ze_na()
  
  
  ## 4. Maternal (1x) great-grandparents
  mmm <- opop %>% filter(pid == mm) %>% pull(mom) %>% ze_na()
  mmf <- opop %>% filter(pid == mm) %>% pull(pop) %>% ze_na()
  mfm <- opop %>% filter(pid == mf) %>% pull(mom) %>% ze_na()
  mff <- opop %>% filter(pid == mf) %>% pull(pop) %>% ze_na()
  
  ## 4. Paternal (1x) great-grandparents
  fmm <- opop %>% filter(pid == fm) %>% pull(mom) %>% ze_na()
  fmf <- opop %>% filter(pid == fm) %>% pull(pop) %>% ze_na()
  ffm <- opop %>% filter(pid == ff) %>% pull(mom) %>% ze_na() 
  fff <- opop %>% filter(pid == ff) %>% pull(pop) %>% ze_na()
  
  
  ## 5. Maternal (2x) great-great-grandparents
  mmmm <- opop %>% filter(pid == mmm) %>% pull(mom) %>% ze_na()
  mmmf <- opop %>% filter(pid == mmm) %>% pull(pop) %>% ze_na()
  mmfm <- opop %>% filter(pid == mmf) %>% pull(mom) %>% ze_na()
  mmff <- opop %>% filter(pid == mmf) %>% pull(pop) %>% ze_na()
  mfmm <- opop %>% filter(pid == mfm) %>% pull(mom) %>% ze_na()
  mfmf <- opop %>% filter(pid == mfm) %>% pull(pop) %>% ze_na()
  mffm <- opop %>% filter(pid == mff) %>% pull(mom) %>% ze_na()
  mfff <- opop %>% filter(pid == mff) %>% pull(pop) %>% ze_na()
  
  ## 5. Paternal (2x) great-great-grandparents
  fmmm <- opop %>% filter(pid == fmm) %>% pull(mom) %>% ze_na()
  fmmf <- opop %>% filter(pid == fmm) %>% pull(pop) %>% ze_na()
  fmfm <- opop %>% filter(pid == fmf) %>% pull(mom) %>% ze_na()
  fmff <- opop %>% filter(pid == fmf) %>% pull(pop) %>% ze_na()
  ffmm <- opop %>% filter(pid == ffm) %>% pull(mom) %>% ze_na()
  ffmf <- opop %>% filter(pid == ffm) %>% pull(pop) %>% ze_na()
  fffm <- opop %>% filter(pid == fff) %>% pull(mom) %>% ze_na()
  ffff <- opop %>% filter(pid == fff) %>% pull(pop) %>% ze_na()
  
  
  ## 6. Maternal (3x) great-great-great-grandparents
  mmmmm <- opop %>% filter(pid == mmmm) %>% pull(mom) %>% ze_na()
  mmmmf <- opop %>% filter(pid == mmmm) %>% pull(pop) %>% ze_na()
  mmmfm <- opop %>% filter(pid == mmmf) %>% pull(mom) %>% ze_na()
  mmmff <- opop %>% filter(pid == mmmf) %>% pull(pop) %>% ze_na()
  mmfmm <- opop %>% filter(pid == mmfm) %>% pull(mom) %>% ze_na()
  mmfmf <- opop %>% filter(pid == mmfm) %>% pull(pop) %>% ze_na()
  mmffm <- opop %>% filter(pid == mmff) %>% pull(mom) %>% ze_na()
  mmfff <- opop %>% filter(pid == mmff) %>% pull(pop) %>% ze_na()
  mfmmm <- opop %>% filter(pid == mfmm) %>% pull(mom) %>% ze_na()
  mfmmf <- opop %>% filter(pid == mfmm) %>% pull(pop) %>% ze_na()
  mfmfm <- opop %>% filter(pid == mfmf) %>% pull(mom) %>% ze_na()
  mfmff <- opop %>% filter(pid == mfmf) %>% pull(pop) %>% ze_na()
  mffmm <- opop %>% filter(pid == mffm) %>% pull(mom) %>% ze_na()
  mffmf <- opop %>% filter(pid == mffm) %>% pull(pop) %>% ze_na()
  mfffm <- opop %>% filter(pid == mfff) %>% pull(mom) %>% ze_na()
  mffff <- opop %>% filter(pid == mfff) %>% pull(pop) %>% ze_na()
  
  ## 6. Paternal (3x) great-great-great-grandparents
  fmmmm <- opop %>% filter(pid == fmmm) %>% pull(mom) %>% ze_na()
  fmmmf <- opop %>% filter(pid == fmmm) %>% pull(pop) %>% ze_na()
  fmmfm <- opop %>% filter(pid == fmmf) %>% pull(mom) %>% ze_na()
  fmmff <- opop %>% filter(pid == fmmf) %>% pull(pop) %>% ze_na()
  fmfmm <- opop %>% filter(pid == fmfm) %>% pull(mom) %>% ze_na()
  fmfmf <- opop %>% filter(pid == fmfm) %>% pull(pop) %>% ze_na()
  fmffm <- opop %>% filter(pid == fmff) %>% pull(mom) %>% ze_na()
  fmfff <- opop %>% filter(pid == fmff) %>% pull(pop) %>% ze_na()
  ffmmm <- opop %>% filter(pid == ffmm) %>% pull(mom) %>% ze_na()
  ffmmf <- opop %>% filter(pid == ffmm) %>% pull(pop) %>% ze_na()
  ffmfm <- opop %>% filter(pid == ffmf) %>% pull(mom) %>% ze_na()
  ffmff <- opop %>% filter(pid == ffmf) %>% pull(pop) %>% ze_na()
  fffmm <- opop %>% filter(pid == fffm) %>% pull(mom) %>% ze_na()
  fffmf <- opop %>% filter(pid == fffm) %>% pull(pop) %>% ze_na()
  ffffm <- opop %>% filter(pid == ffff) %>% pull(mom) %>% ze_na()
  fffff <- opop %>% filter(pid == ffff) %>% pull(pop) %>% ze_na()
  
  
  ## 7. Maternal (4x) great-great-great-great-grandparents
  mmmmmm <- opop %>% filter(pid == mmmmm) %>% pull(mom) %>% ze_na()
  mmmmmf <- opop %>% filter(pid == mmmmm) %>% pull(pop) %>% ze_na()
  mmmmfm <- opop %>% filter(pid == mmmmf) %>% pull(mom) %>% ze_na()
  mmmmff <- opop %>% filter(pid == mmmmf) %>% pull(pop) %>% ze_na()
  mmmfmm <- opop %>% filter(pid == mmmfm) %>% pull(mom) %>% ze_na()
  mmmfmf <- opop %>% filter(pid == mmmfm) %>% pull(pop) %>% ze_na()
  mmmffm <- opop %>% filter(pid == mmmff) %>% pull(mom) %>% ze_na()
  mmmfff <- opop %>% filter(pid == mmmff) %>% pull(pop) %>% ze_na()
  mmfmmm <- opop %>% filter(pid == mmfmm) %>% pull(mom) %>% ze_na()
  mmfmmf <- opop %>% filter(pid == mmfmm) %>% pull(pop) %>% ze_na()
  mmfmfm <- opop %>% filter(pid == mmfmf) %>% pull(mom) %>% ze_na()
  mmfmff <- opop %>% filter(pid == mmfmf) %>% pull(pop) %>% ze_na()
  mmffmm <- opop %>% filter(pid == mmffm) %>% pull(mom) %>% ze_na()
  mmffmf <- opop %>% filter(pid == mmffm) %>% pull(pop) %>% ze_na()
  mmfffm <- opop %>% filter(pid == mmfff) %>% pull(mom) %>% ze_na()
  mmffff <- opop %>% filter(pid == mmfff) %>% pull(pop) %>% ze_na()
  mfmmmm <- opop %>% filter(pid == mfmmm) %>% pull(mom) %>% ze_na()
  mfmmmf <- opop %>% filter(pid == mfmmm) %>% pull(pop) %>% ze_na()
  mfmmfm <- opop %>% filter(pid == mfmmf) %>% pull(mom) %>% ze_na()
  mfmmff <- opop %>% filter(pid == mfmmf) %>% pull(pop) %>% ze_na()
  mfmfmm <- opop %>% filter(pid == mfmfm) %>% pull(mom) %>% ze_na()
  mfmfmf <- opop %>% filter(pid == mfmfm) %>% pull(pop) %>% ze_na()
  mfmffm <- opop %>% filter(pid == mfmff) %>% pull(mom) %>% ze_na()
  mfmfff <- opop %>% filter(pid == mfmff) %>% pull(pop) %>% ze_na()
  mffmmm <- opop %>% filter(pid == mffmm) %>% pull(mom) %>% ze_na()
  mffmmf <- opop %>% filter(pid == mffmm) %>% pull(pop) %>% ze_na()
  mffmfm <- opop %>% filter(pid == mffmf) %>% pull(mom) %>% ze_na()
  mffmff <- opop %>% filter(pid == mffmf) %>% pull(pop) %>% ze_na()
  mfffmm <- opop %>% filter(pid == mfffm) %>% pull(mom) %>% ze_na()
  mfffmf <- opop %>% filter(pid == mfffm) %>% pull(pop) %>% ze_na()
  mffffm <- opop %>% filter(pid == mffff) %>% pull(mom) %>% ze_na()
  mfffff <- opop %>% filter(pid == mffff) %>% pull(pop) %>% ze_na()
  
  
  ## 7. Paternal (4x) great-great-great-great-grandparents
  fmmmmm <- opop %>% filter(pid == fmmmm) %>% pull(mom) %>% ze_na()
  fmmmmf <- opop %>% filter(pid == fmmmm) %>% pull(pop) %>% ze_na()
  fmmmfm <- opop %>% filter(pid == fmmmf) %>% pull(mom) %>% ze_na()
  fmmmff <- opop %>% filter(pid == fmmmf) %>% pull(pop) %>% ze_na()
  fmmfmm <- opop %>% filter(pid == fmmfm) %>% pull(mom) %>% ze_na()
  fmmfmf <- opop %>% filter(pid == fmmfm) %>% pull(pop) %>% ze_na()
  fmmffm <- opop %>% filter(pid == fmmff) %>% pull(mom) %>% ze_na()
  fmmfff <- opop %>% filter(pid == fmmff) %>% pull(pop) %>% ze_na()
  fmfmmm <- opop %>% filter(pid == fmfmm) %>% pull(mom) %>% ze_na()
  fmfmmf <- opop %>% filter(pid == fmfmm) %>% pull(pop) %>% ze_na()
  fmfmfm <- opop %>% filter(pid == fmfmf) %>% pull(mom) %>% ze_na()
  fmfmff <- opop %>% filter(pid == fmfmf) %>% pull(pop) %>% ze_na()
  fmffmm <- opop %>% filter(pid == fmffm) %>% pull(mom) %>% ze_na()
  fmffmf <- opop %>% filter(pid == fmffm) %>% pull(pop) %>% ze_na()
  fmfffm <- opop %>% filter(pid == fmfff) %>% pull(mom) %>% ze_na()
  fmffff <- opop %>% filter(pid == fmfff) %>% pull(pop) %>% ze_na()
  ffmmmm <- opop %>% filter(pid == ffmmm) %>% pull(mom) %>% ze_na()
  ffmmmf <- opop %>% filter(pid == ffmmm) %>% pull(pop) %>% ze_na()
  ffmmfm <- opop %>% filter(pid == ffmmf) %>% pull(mom) %>% ze_na()
  ffmmff <- opop %>% filter(pid == ffmmf) %>% pull(pop) %>% ze_na()
  ffmfmm <- opop %>% filter(pid == ffmfm) %>% pull(mom) %>% ze_na()
  ffmfmf <- opop %>% filter(pid == ffmfm) %>% pull(pop) %>% ze_na()
  ffmffm <- opop %>% filter(pid == ffmff) %>% pull(mom) %>% ze_na()
  ffmfff <- opop %>% filter(pid == ffmff) %>% pull(pop) %>% ze_na()
  fffmmm <- opop %>% filter(pid == fffmm) %>% pull(mom) %>% ze_na()
  fffmmf <- opop %>% filter(pid == fffmm) %>% pull(pop) %>% ze_na()
  fffmfm <- opop %>% filter(pid == fffmf) %>% pull(mom) %>% ze_na()
  fffmff <- opop %>% filter(pid == fffmf) %>% pull(pop) %>% ze_na()
  ffffmm <- opop %>% filter(pid == ffffm) %>% pull(mom) %>% ze_na()
  ffffmf <- opop %>% filter(pid == ffffm) %>% pull(pop) %>% ze_na()
  fffffm <- opop %>% filter(pid == fffff) %>% pull(mom) %>% ze_na()
  ffffff <- opop %>% filter(pid == fffff) %>% pull(pop) %>% ze_na()
  
  
  ## 8. Maternal (5x) great-great-great-great-great-grandparents
  mmmmmmm <- opop %>% filter(pid == mmmmmm) %>% pull(mom) %>% ze_na()
  mmmmmmf <- opop %>% filter(pid == mmmmmm) %>% pull(pop) %>% ze_na()
  mmmmmfm <- opop %>% filter(pid == mmmmmf) %>% pull(mom) %>% ze_na()
  mmmmmff <- opop %>% filter(pid == mmmmmf) %>% pull(pop) %>% ze_na()
  mmmmfmm <- opop %>% filter(pid == mmmmfm) %>% pull(mom) %>% ze_na()
  mmmmfmf <- opop %>% filter(pid == mmmmfm) %>% pull(pop) %>% ze_na()
  mmmmffm <- opop %>% filter(pid == mmmmff) %>% pull(mom) %>% ze_na()
  mmmmfff <- opop %>% filter(pid == mmmmff) %>% pull(pop) %>% ze_na()
  mmmfmmm <- opop %>% filter(pid == mmmfmm) %>% pull(mom) %>% ze_na()
  mmmfmmf <- opop %>% filter(pid == mmmfmm) %>% pull(pop) %>% ze_na()
  mmmfmfm <- opop %>% filter(pid == mmmfmf) %>% pull(mom) %>% ze_na()
  mmmfmff <- opop %>% filter(pid == mmmfmf) %>% pull(pop) %>% ze_na()
  mmmffmm <- opop %>% filter(pid == mmmffm) %>% pull(mom) %>% ze_na()
  mmmffmf <- opop %>% filter(pid == mmmffm) %>% pull(pop) %>% ze_na()
  mmmfffm <- opop %>% filter(pid == mmmfff) %>% pull(mom) %>% ze_na()
  mmmffff <- opop %>% filter(pid == mmmfff) %>% pull(pop) %>% ze_na()
  mmfmmmm <- opop %>% filter(pid == mmfmmm) %>% pull(mom) %>% ze_na()
  mmfmmmf <- opop %>% filter(pid == mmfmmm) %>% pull(pop) %>% ze_na()
  mmfmmfm <- opop %>% filter(pid == mmfmmf) %>% pull(mom) %>% ze_na()
  mmfmmff <- opop %>% filter(pid == mmfmmf) %>% pull(pop) %>% ze_na()
  mmfmfmm <- opop %>% filter(pid == mmfmfm) %>% pull(mom) %>% ze_na()
  mmfmfmf <- opop %>% filter(pid == mmfmfm) %>% pull(pop) %>% ze_na()
  mmfmffm <- opop %>% filter(pid == mmfmff) %>% pull(mom) %>% ze_na()
  mmfmfff <- opop %>% filter(pid == mmfmff) %>% pull(pop) %>% ze_na()
  mmffmmm <- opop %>% filter(pid == mmffmm) %>% pull(mom) %>% ze_na()
  mmffmmf <- opop %>% filter(pid == mmffmm) %>% pull(pop) %>% ze_na()
  mmffmfm <- opop %>% filter(pid == mmffmf) %>% pull(mom) %>% ze_na()
  mmffmff <- opop %>% filter(pid == mmffmf) %>% pull(pop) %>% ze_na()
  mmfffmm <- opop %>% filter(pid == mmfffm) %>% pull(mom) %>% ze_na()
  mmfffmf <- opop %>% filter(pid == mmfffm) %>% pull(pop) %>% ze_na()
  mmffffm <- opop %>% filter(pid == mmffff) %>% pull(mom) %>% ze_na()
  mmfffff <- opop %>% filter(pid == mmffff) %>% pull(pop) %>% ze_na()
  mfmmmmm <- opop %>% filter(pid == mfmmmm) %>% pull(mom) %>% ze_na()
  mfmmmmf <- opop %>% filter(pid == mfmmmm) %>% pull(pop) %>% ze_na()
  mfmmmfm <- opop %>% filter(pid == mfmmmf) %>% pull(mom) %>% ze_na()
  mfmmmff <- opop %>% filter(pid == mfmmmf) %>% pull(pop) %>% ze_na()
  mfmmfmm <- opop %>% filter(pid == mfmmfm) %>% pull(mom) %>% ze_na()
  mfmmfmf <- opop %>% filter(pid == mfmmfm) %>% pull(pop) %>% ze_na()
  mfmmffm <- opop %>% filter(pid == mfmmff) %>% pull(mom) %>% ze_na()
  mfmmfff <- opop %>% filter(pid == mfmmff) %>% pull(pop) %>% ze_na()
  mfmfmmm <- opop %>% filter(pid == mfmfmm) %>% pull(mom) %>% ze_na()
  mfmfmmf <- opop %>% filter(pid == mfmfmm) %>% pull(pop) %>% ze_na()
  mfmfmfm <- opop %>% filter(pid == mfmfmf) %>% pull(mom) %>% ze_na()
  mfmfmff <- opop %>% filter(pid == mfmfmf) %>% pull(pop) %>% ze_na()
  mfmffmm <- opop %>% filter(pid == mfmffm) %>% pull(mom) %>% ze_na()
  mfmffmf <- opop %>% filter(pid == mfmffm) %>% pull(pop) %>% ze_na()
  mfmfffm <- opop %>% filter(pid == mfmfff) %>% pull(mom) %>% ze_na()
  mfmffff <- opop %>% filter(pid == mfmfff) %>% pull(pop) %>% ze_na()
  mffmmmm <- opop %>% filter(pid == mffmmm) %>% pull(mom) %>% ze_na()
  mffmmmf <- opop %>% filter(pid == mffmmm) %>% pull(pop) %>% ze_na()
  mffmmfm <- opop %>% filter(pid == mffmmf) %>% pull(mom) %>% ze_na()
  mffmmff <- opop %>% filter(pid == mffmmf) %>% pull(pop) %>% ze_na()
  mffmfmm <- opop %>% filter(pid == mffmfm) %>% pull(mom) %>% ze_na()
  mffmfmf <- opop %>% filter(pid == mffmfm) %>% pull(pop) %>% ze_na()
  mffmffm <- opop %>% filter(pid == mffmff) %>% pull(mom) %>% ze_na()
  mffmfff <- opop %>% filter(pid == mffmff) %>% pull(pop) %>% ze_na()
  mfffmmm <- opop %>% filter(pid == mfffmm) %>% pull(mom) %>% ze_na()
  mfffmmf <- opop %>% filter(pid == mfffmm) %>% pull(pop) %>% ze_na()
  mfffmfm <- opop %>% filter(pid == mfffmf) %>% pull(mom) %>% ze_na()
  mfffmff <- opop %>% filter(pid == mfffmf) %>% pull(pop) %>% ze_na()
  mffffmm <- opop %>% filter(pid == mffffm) %>% pull(mom) %>% ze_na()
  mffffmf <- opop %>% filter(pid == mffffm) %>% pull(pop) %>% ze_na()
  mfffffm <- opop %>% filter(pid == mfffff) %>% pull(mom) %>% ze_na()
  mffffff <- opop %>% filter(pid == mfffff) %>% pull(pop) %>% ze_na()
  
  
  ## 8. Paternal (5x) great-great-great-great-great-grandparents
  fmmmmmm <- opop %>% filter(pid == fmmmmm) %>% pull(mom) %>% ze_na()
  fmmmmmf <- opop %>% filter(pid == fmmmmm) %>% pull(pop) %>% ze_na()
  fmmmmfm <- opop %>% filter(pid == fmmmmf) %>% pull(mom) %>% ze_na()
  fmmmmff <- opop %>% filter(pid == fmmmmf) %>% pull(pop) %>% ze_na()
  fmmmfmm <- opop %>% filter(pid == fmmmfm) %>% pull(mom) %>% ze_na()
  fmmmfmf <- opop %>% filter(pid == fmmmfm) %>% pull(pop) %>% ze_na()
  fmmmffm <- opop %>% filter(pid == fmmmff) %>% pull(mom) %>% ze_na()
  fmmmfff <- opop %>% filter(pid == fmmmff) %>% pull(pop) %>% ze_na()
  fmmfmmm <- opop %>% filter(pid == fmmfmm) %>% pull(mom) %>% ze_na()
  fmmfmmf <- opop %>% filter(pid == fmmfmm) %>% pull(pop) %>% ze_na()
  fmmfmfm <- opop %>% filter(pid == fmmfmf) %>% pull(mom) %>% ze_na()
  fmmfmff <- opop %>% filter(pid == fmmfmf) %>% pull(pop) %>% ze_na()
  fmmffmm <- opop %>% filter(pid == fmmffm) %>% pull(mom) %>% ze_na()
  fmmffmf <- opop %>% filter(pid == fmmffm) %>% pull(pop) %>% ze_na()
  fmmfffm <- opop %>% filter(pid == fmmfff) %>% pull(mom) %>% ze_na()
  fmmffff <- opop %>% filter(pid == fmmfff) %>% pull(pop) %>% ze_na()
  fmfmmmm <- opop %>% filter(pid == fmfmmm) %>% pull(mom) %>% ze_na()
  fmfmmmf <- opop %>% filter(pid == fmfmmm) %>% pull(pop) %>% ze_na()
  fmfmmfm <- opop %>% filter(pid == fmfmmf) %>% pull(mom) %>% ze_na()
  fmfmmff <- opop %>% filter(pid == fmfmmf) %>% pull(pop) %>% ze_na()
  fmfmfmm <- opop %>% filter(pid == fmfmfm) %>% pull(mom) %>% ze_na()
  fmfmfmf <- opop %>% filter(pid == fmfmfm) %>% pull(pop) %>% ze_na()
  fmfmffm <- opop %>% filter(pid == fmfmff) %>% pull(mom) %>% ze_na()
  fmfmfff <- opop %>% filter(pid == fmfmff) %>% pull(pop) %>% ze_na()
  fmffmmm <- opop %>% filter(pid == fmffmm) %>% pull(mom) %>% ze_na()
  fmffmmf <- opop %>% filter(pid == fmffmm) %>% pull(pop) %>% ze_na()
  fmffmfm <- opop %>% filter(pid == fmffmf) %>% pull(mom) %>% ze_na()
  fmffmff <- opop %>% filter(pid == fmffmf) %>% pull(pop) %>% ze_na()
  fmfffmm <- opop %>% filter(pid == fmfffm) %>% pull(mom) %>% ze_na()
  fmfffmf <- opop %>% filter(pid == fmfffm) %>% pull(pop) %>% ze_na()
  fmffffm <- opop %>% filter(pid == fmffff) %>% pull(mom) %>% ze_na()
  fmfffff <- opop %>% filter(pid == fmffff) %>% pull(pop) %>% ze_na()
  ffmmmmm <- opop %>% filter(pid == ffmmmm) %>% pull(mom) %>% ze_na()
  ffmmmmf <- opop %>% filter(pid == ffmmmm) %>% pull(pop) %>% ze_na()
  ffmmmfm <- opop %>% filter(pid == ffmmmf) %>% pull(mom) %>% ze_na()
  ffmmmff <- opop %>% filter(pid == ffmmmf) %>% pull(pop) %>% ze_na()
  ffmmfmm <- opop %>% filter(pid == ffmmfm) %>% pull(mom) %>% ze_na()
  ffmmfmf <- opop %>% filter(pid == ffmmfm) %>% pull(pop) %>% ze_na()
  ffmmffm <- opop %>% filter(pid == ffmmff) %>% pull(mom) %>% ze_na()
  ffmmfff <- opop %>% filter(pid == ffmmff) %>% pull(pop) %>% ze_na()
  ffmfmmm <- opop %>% filter(pid == ffmfmm) %>% pull(mom) %>% ze_na()
  ffmfmmf <- opop %>% filter(pid == ffmfmm) %>% pull(pop) %>% ze_na()
  ffmfmfm <- opop %>% filter(pid == ffmfmf) %>% pull(mom) %>% ze_na()
  ffmfmff <- opop %>% filter(pid == ffmfmf) %>% pull(pop) %>% ze_na()
  ffmffmm <- opop %>% filter(pid == ffmffm) %>% pull(mom) %>% ze_na()
  ffmffmf <- opop %>% filter(pid == ffmffm) %>% pull(pop) %>% ze_na()
  ffmfffm <- opop %>% filter(pid == ffmfff) %>% pull(mom) %>% ze_na()
  ffmffff <- opop %>% filter(pid == ffmfff) %>% pull(pop) %>% ze_na()
  fffmmmm <- opop %>% filter(pid == fffmmm) %>% pull(mom) %>% ze_na()
  fffmmmf <- opop %>% filter(pid == fffmmm) %>% pull(pop) %>% ze_na()
  fffmmfm <- opop %>% filter(pid == fffmmf) %>% pull(mom) %>% ze_na()
  fffmmff <- opop %>% filter(pid == fffmmf) %>% pull(pop) %>% ze_na()
  fffmfmm <- opop %>% filter(pid == fffmfm) %>% pull(mom) %>% ze_na()
  fffmfmf <- opop %>% filter(pid == fffmfm) %>% pull(pop) %>% ze_na()
  fffmffm <- opop %>% filter(pid == fffmff) %>% pull(mom) %>% ze_na()
  fffmfff <- opop %>% filter(pid == fffmff) %>% pull(pop) %>% ze_na()
  ffffmmm <- opop %>% filter(pid == ffffmm) %>% pull(mom) %>% ze_na()
  ffffmmf <- opop %>% filter(pid == ffffmm) %>% pull(pop) %>% ze_na()
  ffffmfm <- opop %>% filter(pid == ffffmf) %>% pull(mom) %>% ze_na()
  ffffmff <- opop %>% filter(pid == ffffmf) %>% pull(pop) %>% ze_na()
  fffffmm <- opop %>% filter(pid == fffffm) %>% pull(mom) %>% ze_na()
  fffffmf <- opop %>% filter(pid == fffffm) %>% pull(pop) %>% ze_na()
  ffffffm <- opop %>% filter(pid == ffffff) %>% pull(mom) %>% ze_na()
  fffffff <- opop %>% filter(pid == ffffff) %>% pull(pop) %>% ze_na()
  
  
  ## 8. Maternal (6x) great-great-great-great-great-great-grandparents
  mmmmmmmm <- opop %>% filter(pid == mmmmmmm) %>% pull(mom) %>% ze_na()
  mmmmmmmf <- opop %>% filter(pid == mmmmmmm) %>% pull(pop) %>% ze_na()
  mmmmmmfm <- opop %>% filter(pid == mmmmmmf) %>% pull(mom) %>% ze_na()
  mmmmmmff <- opop %>% filter(pid == mmmmmmf) %>% pull(pop) %>% ze_na()
  mmmmmfmm <- opop %>% filter(pid == mmmmmfm) %>% pull(mom) %>% ze_na()
  mmmmmfmf <- opop %>% filter(pid == mmmmmfm) %>% pull(pop) %>% ze_na()
  mmmmmffm <- opop %>% filter(pid == mmmmmff) %>% pull(mom) %>% ze_na()
  mmmmmfff <- opop %>% filter(pid == mmmmmff) %>% pull(pop) %>% ze_na()
  mmmmfmmm <- opop %>% filter(pid == mmmmfmm) %>% pull(mom) %>% ze_na()
  mmmmfmmf <- opop %>% filter(pid == mmmmfmm) %>% pull(pop) %>% ze_na()
  mmmmfmfm <- opop %>% filter(pid == mmmmfmf) %>% pull(mom) %>% ze_na()
  mmmmfmff <- opop %>% filter(pid == mmmmfmf) %>% pull(pop) %>% ze_na()
  mmmmffmm <- opop %>% filter(pid == mmmmffm) %>% pull(mom) %>% ze_na()
  mmmmffmf <- opop %>% filter(pid == mmmmffm) %>% pull(pop) %>% ze_na()
  mmmmfffm <- opop %>% filter(pid == mmmmfff) %>% pull(mom) %>% ze_na()
  mmmmffff <- opop %>% filter(pid == mmmmfff) %>% pull(pop) %>% ze_na()
  mmmfmmmm <- opop %>% filter(pid == mmmfmmm) %>% pull(mom) %>% ze_na()
  mmmfmmmf <- opop %>% filter(pid == mmmfmmm) %>% pull(pop) %>% ze_na()
  mmmfmmfm <- opop %>% filter(pid == mmmfmmf) %>% pull(mom) %>% ze_na()
  mmmfmmff <- opop %>% filter(pid == mmmfmmf) %>% pull(pop) %>% ze_na()
  mmmfmfmm <- opop %>% filter(pid == mmmfmfm) %>% pull(mom) %>% ze_na()
  mmmfmfmf <- opop %>% filter(pid == mmmfmfm) %>% pull(pop) %>% ze_na()
  mmmfmffm <- opop %>% filter(pid == mmmfmff) %>% pull(mom) %>% ze_na()
  mmmfmfff <- opop %>% filter(pid == mmmfmff) %>% pull(pop) %>% ze_na()
  mmmffmmm <- opop %>% filter(pid == mmmffmm) %>% pull(mom) %>% ze_na()
  mmmffmmf <- opop %>% filter(pid == mmmffmm) %>% pull(pop) %>% ze_na()
  mmmffmfm <- opop %>% filter(pid == mmmffmf) %>% pull(mom) %>% ze_na()
  mmmffmff <- opop %>% filter(pid == mmmffmf) %>% pull(pop) %>% ze_na()
  mmmfffmm <- opop %>% filter(pid == mmmfffm) %>% pull(mom) %>% ze_na()
  mmmfffmf <- opop %>% filter(pid == mmmfffm) %>% pull(pop) %>% ze_na()
  mmmffffm <- opop %>% filter(pid == mmmffff) %>% pull(mom) %>% ze_na()
  mmmfffff <- opop %>% filter(pid == mmmffff) %>% pull(pop) %>% ze_na()
  mmfmmmmm <- opop %>% filter(pid == mmfmmmm) %>% pull(mom) %>% ze_na()
  mmfmmmmf <- opop %>% filter(pid == mmfmmmm) %>% pull(pop) %>% ze_na()
  mmfmmmfm <- opop %>% filter(pid == mmfmmmf) %>% pull(mom) %>% ze_na()
  mmfmmmff <- opop %>% filter(pid == mmfmmmf) %>% pull(pop) %>% ze_na()
  mmfmmfmm <- opop %>% filter(pid == mmfmmfm) %>% pull(mom) %>% ze_na()
  mmfmmfmf <- opop %>% filter(pid == mmfmmfm) %>% pull(pop) %>% ze_na()
  mmfmmffm <- opop %>% filter(pid == mmfmmff) %>% pull(mom) %>% ze_na()
  mmfmmfff <- opop %>% filter(pid == mmfmmff) %>% pull(pop) %>% ze_na()
  mmfmfmmm <- opop %>% filter(pid == mmfmfmm) %>% pull(mom) %>% ze_na()
  mmfmfmmf <- opop %>% filter(pid == mmfmfmm) %>% pull(pop) %>% ze_na()
  mmfmfmfm <- opop %>% filter(pid == mmfmfmf) %>% pull(mom) %>% ze_na()
  mmfmfmff <- opop %>% filter(pid == mmfmfmf) %>% pull(pop) %>% ze_na()
  mmfmffmm <- opop %>% filter(pid == mmfmffm) %>% pull(mom) %>% ze_na()
  mmfmffmf <- opop %>% filter(pid == mmfmffm) %>% pull(pop) %>% ze_na()
  mmfmfffm <- opop %>% filter(pid == mmfmfff) %>% pull(mom) %>% ze_na()
  mmfmffff <- opop %>% filter(pid == mmfmfff) %>% pull(pop) %>% ze_na()
  mmffmmmm <- opop %>% filter(pid == mmffmmm) %>% pull(mom) %>% ze_na()
  mmffmmmf <- opop %>% filter(pid == mmffmmm) %>% pull(pop) %>% ze_na()
  mmffmmfm <- opop %>% filter(pid == mmffmmf) %>% pull(mom) %>% ze_na()
  mmffmmff <- opop %>% filter(pid == mmffmmf) %>% pull(pop) %>% ze_na()
  mmffmfmm <- opop %>% filter(pid == mmffmfm) %>% pull(mom) %>% ze_na()
  mmffmfmf <- opop %>% filter(pid == mmffmfm) %>% pull(pop) %>% ze_na()
  mmffmffm <- opop %>% filter(pid == mmffmff) %>% pull(mom) %>% ze_na()
  mmffmfff <- opop %>% filter(pid == mmffmff) %>% pull(pop) %>% ze_na()
  mmfffmmm <- opop %>% filter(pid == mmfffmm) %>% pull(mom) %>% ze_na()
  mmfffmmf <- opop %>% filter(pid == mmfffmm) %>% pull(pop) %>% ze_na()
  mmfffmfm <- opop %>% filter(pid == mmfffmf) %>% pull(mom) %>% ze_na()
  mmfffmff <- opop %>% filter(pid == mmfffmf) %>% pull(pop) %>% ze_na()
  mmffffmm <- opop %>% filter(pid == mmffffm) %>% pull(mom) %>% ze_na()
  mmffffmf <- opop %>% filter(pid == mmffffm) %>% pull(pop) %>% ze_na()
  mmfffffm <- opop %>% filter(pid == mmfffff) %>% pull(mom) %>% ze_na()
  mmffffff <- opop %>% filter(pid == mmfffff) %>% pull(pop) %>% ze_na()
  mfmmmmmm <- opop %>% filter(pid == mfmmmmm) %>% pull(mom) %>% ze_na()
  mfmmmmmf <- opop %>% filter(pid == mfmmmmm) %>% pull(pop) %>% ze_na()
  mfmmmmfm <- opop %>% filter(pid == mfmmmmf) %>% pull(mom) %>% ze_na()
  mfmmmmff <- opop %>% filter(pid == mfmmmmf) %>% pull(pop) %>% ze_na()
  mfmmmfmm <- opop %>% filter(pid == mfmmmfm) %>% pull(mom) %>% ze_na()
  mfmmmfmf <- opop %>% filter(pid == mfmmmfm) %>% pull(pop) %>% ze_na()
  mfmmmffm <- opop %>% filter(pid == mfmmmff) %>% pull(mom) %>% ze_na()
  mfmmmfff <- opop %>% filter(pid == mfmmmff) %>% pull(pop) %>% ze_na()
  mfmmfmmm <- opop %>% filter(pid == mfmmfmm) %>% pull(mom) %>% ze_na()
  mfmmfmmf <- opop %>% filter(pid == mfmmfmm) %>% pull(pop) %>% ze_na()
  mfmmfmfm <- opop %>% filter(pid == mfmmfmf) %>% pull(mom) %>% ze_na()
  mfmmfmff <- opop %>% filter(pid == mfmmfmf) %>% pull(pop) %>% ze_na()
  mfmmffmm <- opop %>% filter(pid == mfmmffm) %>% pull(mom) %>% ze_na()
  mfmmffmf <- opop %>% filter(pid == mfmmffm) %>% pull(pop) %>% ze_na()
  mfmmfffm <- opop %>% filter(pid == mfmmfff) %>% pull(mom) %>% ze_na()
  mfmmffff <- opop %>% filter(pid == mfmmfff) %>% pull(pop) %>% ze_na()
  mfmfmmmm <- opop %>% filter(pid == mfmfmmm) %>% pull(mom) %>% ze_na()
  mfmfmmmf <- opop %>% filter(pid == mfmfmmm) %>% pull(pop) %>% ze_na()
  mfmfmmfm <- opop %>% filter(pid == mfmfmmf) %>% pull(mom) %>% ze_na()
  mfmfmmff <- opop %>% filter(pid == mfmfmmf) %>% pull(pop) %>% ze_na()
  mfmfmfmm <- opop %>% filter(pid == mfmfmfm) %>% pull(mom) %>% ze_na()
  mfmfmfmf <- opop %>% filter(pid == mfmfmfm) %>% pull(pop) %>% ze_na()
  mfmfmffm <- opop %>% filter(pid == mfmfmff) %>% pull(mom) %>% ze_na()
  mfmfmfff <- opop %>% filter(pid == mfmfmff) %>% pull(pop) %>% ze_na()
  mfmffmmm <- opop %>% filter(pid == mfmffmm) %>% pull(mom) %>% ze_na()
  mfmffmmf <- opop %>% filter(pid == mfmffmm) %>% pull(pop) %>% ze_na()
  mfmffmfm <- opop %>% filter(pid == mfmffmf) %>% pull(mom) %>% ze_na()
  mfmffmff <- opop %>% filter(pid == mfmffmf) %>% pull(pop) %>% ze_na()
  mfmfffmm <- opop %>% filter(pid == mfmfffm) %>% pull(mom) %>% ze_na()
  mfmfffmf <- opop %>% filter(pid == mfmfffm) %>% pull(pop) %>% ze_na()
  mfmffffm <- opop %>% filter(pid == mfmffff) %>% pull(mom) %>% ze_na()
  mfmfffff <- opop %>% filter(pid == mfmffff) %>% pull(pop) %>% ze_na()
  mffmmmmm <- opop %>% filter(pid == mffmmmm) %>% pull(mom) %>% ze_na()
  mffmmmmf <- opop %>% filter(pid == mffmmmm) %>% pull(pop) %>% ze_na()
  mffmmmfm <- opop %>% filter(pid == mffmmmf) %>% pull(mom) %>% ze_na()
  mffmmmff <- opop %>% filter(pid == mffmmmf) %>% pull(pop) %>% ze_na()
  mffmmfmm <- opop %>% filter(pid == mffmmfm) %>% pull(mom) %>% ze_na()
  mffmmfmf <- opop %>% filter(pid == mffmmfm) %>% pull(pop) %>% ze_na()
  mffmmffm <- opop %>% filter(pid == mffmmff) %>% pull(mom) %>% ze_na()
  mffmmfff <- opop %>% filter(pid == mffmmff) %>% pull(pop) %>% ze_na()
  mffmfmmm <- opop %>% filter(pid == mffmfmm) %>% pull(mom) %>% ze_na()
  mffmfmmf <- opop %>% filter(pid == mffmfmm) %>% pull(pop) %>% ze_na()
  mffmfmfm <- opop %>% filter(pid == mffmfmf) %>% pull(mom) %>% ze_na()
  mffmfmff <- opop %>% filter(pid == mffmfmf) %>% pull(pop) %>% ze_na()
  mffmffmm <- opop %>% filter(pid == mffmffm) %>% pull(mom) %>% ze_na()
  mffmffmf <- opop %>% filter(pid == mffmffm) %>% pull(pop) %>% ze_na()
  mffmfffm <- opop %>% filter(pid == mffmfff) %>% pull(mom) %>% ze_na()
  mffmffff <- opop %>% filter(pid == mffmfff) %>% pull(pop) %>% ze_na()
  mfffmmmm <- opop %>% filter(pid == mfffmmm) %>% pull(mom) %>% ze_na()
  mfffmmmf <- opop %>% filter(pid == mfffmmm) %>% pull(pop) %>% ze_na()
  mfffmmfm <- opop %>% filter(pid == mfffmmf) %>% pull(mom) %>% ze_na()
  mfffmmff <- opop %>% filter(pid == mfffmmf) %>% pull(pop) %>% ze_na()
  mfffmfmm <- opop %>% filter(pid == mfffmfm) %>% pull(mom) %>% ze_na()
  mfffmfmf <- opop %>% filter(pid == mfffmfm) %>% pull(pop) %>% ze_na()
  mfffmffm <- opop %>% filter(pid == mfffmff) %>% pull(mom) %>% ze_na()
  mfffmfff <- opop %>% filter(pid == mfffmff) %>% pull(pop) %>% ze_na()
  mffffmmm <- opop %>% filter(pid == mffffmm) %>% pull(mom) %>% ze_na()
  mffffmmf <- opop %>% filter(pid == mffffmm) %>% pull(pop) %>% ze_na()
  mffffmfm <- opop %>% filter(pid == mffffmf) %>% pull(mom) %>% ze_na()
  mffffmff <- opop %>% filter(pid == mffffmf) %>% pull(pop) %>% ze_na()
  mfffffmm <- opop %>% filter(pid == mfffffm) %>% pull(mom) %>% ze_na()
  mfffffmf <- opop %>% filter(pid == mfffffm) %>% pull(pop) %>% ze_na()
  mffffffm <- opop %>% filter(pid == mffffff) %>% pull(mom) %>% ze_na()
  mfffffff <- opop %>% filter(pid == mffffff) %>% pull(pop) %>% ze_na()
  
  
  ## 8. Paternal (6x) great-great-great-great-great-great-grandparents
  fmmmmmmm <- opop %>% filter(pid == fmmmmmm) %>% pull(mom) %>% ze_na()
  fmmmmmmf <- opop %>% filter(pid == fmmmmmm) %>% pull(pop) %>% ze_na()
  fmmmmmfm <- opop %>% filter(pid == fmmmmmf) %>% pull(mom) %>% ze_na()
  fmmmmmff <- opop %>% filter(pid == fmmmmmf) %>% pull(pop) %>% ze_na()
  fmmmmfmm <- opop %>% filter(pid == fmmmmfm) %>% pull(mom) %>% ze_na()
  fmmmmfmf <- opop %>% filter(pid == fmmmmfm) %>% pull(pop) %>% ze_na()
  fmmmmffm <- opop %>% filter(pid == fmmmmff) %>% pull(mom) %>% ze_na()
  fmmmmfff <- opop %>% filter(pid == fmmmmff) %>% pull(pop) %>% ze_na()
  fmmmfmmm <- opop %>% filter(pid == fmmmfmm) %>% pull(mom) %>% ze_na()
  fmmmfmmf <- opop %>% filter(pid == fmmmfmm) %>% pull(pop) %>% ze_na()
  fmmmfmfm <- opop %>% filter(pid == fmmmfmf) %>% pull(mom) %>% ze_na()
  fmmmfmff <- opop %>% filter(pid == fmmmfmf) %>% pull(pop) %>% ze_na()
  fmmmffmm <- opop %>% filter(pid == fmmmffm) %>% pull(mom) %>% ze_na()
  fmmmffmf <- opop %>% filter(pid == fmmmffm) %>% pull(pop) %>% ze_na()
  fmmmfffm <- opop %>% filter(pid == fmmmfff) %>% pull(mom) %>% ze_na()
  fmmmffff <- opop %>% filter(pid == fmmmfff) %>% pull(pop) %>% ze_na()
  fmmfmmmm <- opop %>% filter(pid == fmmfmmm) %>% pull(mom) %>% ze_na()
  fmmfmmmf <- opop %>% filter(pid == fmmfmmm) %>% pull(pop) %>% ze_na()
  fmmfmmfm <- opop %>% filter(pid == fmmfmmf) %>% pull(mom) %>% ze_na()
  fmmfmmff <- opop %>% filter(pid == fmmfmmf) %>% pull(pop) %>% ze_na()
  fmmfmfmm <- opop %>% filter(pid == fmmfmfm) %>% pull(mom) %>% ze_na()
  fmmfmfmf <- opop %>% filter(pid == fmmfmfm) %>% pull(pop) %>% ze_na()
  fmmfmffm <- opop %>% filter(pid == fmmfmff) %>% pull(mom) %>% ze_na()
  fmmfmfff <- opop %>% filter(pid == fmmfmff) %>% pull(pop) %>% ze_na()
  fmmffmmm <- opop %>% filter(pid == fmmffmm) %>% pull(mom) %>% ze_na()
  fmmffmmf <- opop %>% filter(pid == fmmffmm) %>% pull(pop) %>% ze_na()
  fmmffmfm <- opop %>% filter(pid == fmmffmf) %>% pull(mom) %>% ze_na()
  fmmffmff <- opop %>% filter(pid == fmmffmf) %>% pull(pop) %>% ze_na()
  fmmfffmm <- opop %>% filter(pid == fmmfffm) %>% pull(mom) %>% ze_na()
  fmmfffmf <- opop %>% filter(pid == fmmfffm) %>% pull(pop) %>% ze_na()
  fmmffffm <- opop %>% filter(pid == fmmffff) %>% pull(mom) %>% ze_na()
  fmmfffff <- opop %>% filter(pid == fmmffff) %>% pull(pop) %>% ze_na()
  fmfmmmmm <- opop %>% filter(pid == fmfmmmm) %>% pull(mom) %>% ze_na()
  fmfmmmmf <- opop %>% filter(pid == fmfmmmm) %>% pull(pop) %>% ze_na()
  fmfmmmfm <- opop %>% filter(pid == fmfmmmf) %>% pull(mom) %>% ze_na()
  fmfmmmff <- opop %>% filter(pid == fmfmmmf) %>% pull(pop) %>% ze_na()
  fmfmmfmm <- opop %>% filter(pid == fmfmmfm) %>% pull(mom) %>% ze_na()
  fmfmmfmf <- opop %>% filter(pid == fmfmmfm) %>% pull(pop) %>% ze_na()
  fmfmmffm <- opop %>% filter(pid == fmfmmff) %>% pull(mom) %>% ze_na()
  fmfmmfff <- opop %>% filter(pid == fmfmmff) %>% pull(pop) %>% ze_na()
  fmfmfmmm <- opop %>% filter(pid == fmfmfmm) %>% pull(mom) %>% ze_na()
  fmfmfmmf <- opop %>% filter(pid == fmfmfmm) %>% pull(pop) %>% ze_na()
  fmfmfmfm <- opop %>% filter(pid == fmfmfmf) %>% pull(mom) %>% ze_na()
  fmfmfmff <- opop %>% filter(pid == fmfmfmf) %>% pull(pop) %>% ze_na()
  fmfmffmm <- opop %>% filter(pid == fmfmffm) %>% pull(mom) %>% ze_na()
  fmfmffmf <- opop %>% filter(pid == fmfmffm) %>% pull(pop) %>% ze_na()
  fmfmfffm <- opop %>% filter(pid == fmfmfff) %>% pull(mom) %>% ze_na()
  fmfmffff <- opop %>% filter(pid == fmfmfff) %>% pull(pop) %>% ze_na()
  fmffmmmm <- opop %>% filter(pid == fmffmmm) %>% pull(mom) %>% ze_na()
  fmffmmmf <- opop %>% filter(pid == fmffmmm) %>% pull(pop) %>% ze_na()
  fmffmmfm <- opop %>% filter(pid == fmffmmf) %>% pull(mom) %>% ze_na()
  fmffmmff <- opop %>% filter(pid == fmffmmf) %>% pull(pop) %>% ze_na()
  fmffmfmm <- opop %>% filter(pid == fmffmfm) %>% pull(mom) %>% ze_na()
  fmffmfmf <- opop %>% filter(pid == fmffmfm) %>% pull(pop) %>% ze_na()
  fmffmffm <- opop %>% filter(pid == fmffmff) %>% pull(mom) %>% ze_na()
  fmffmfff <- opop %>% filter(pid == fmffmff) %>% pull(pop) %>% ze_na()
  fmfffmmm <- opop %>% filter(pid == fmfffmm) %>% pull(mom) %>% ze_na()
  fmfffmmf <- opop %>% filter(pid == fmfffmm) %>% pull(pop) %>% ze_na()
  fmfffmfm <- opop %>% filter(pid == fmfffmf) %>% pull(mom) %>% ze_na()
  fmfffmff <- opop %>% filter(pid == fmfffmf) %>% pull(pop) %>% ze_na()
  fmffffmm <- opop %>% filter(pid == fmffffm) %>% pull(mom) %>% ze_na()
  fmffffmf <- opop %>% filter(pid == fmffffm) %>% pull(pop) %>% ze_na()
  fmfffffm <- opop %>% filter(pid == fmfffff) %>% pull(mom) %>% ze_na()
  fmffffff <- opop %>% filter(pid == fmfffff) %>% pull(pop) %>% ze_na()
  ffmmmmmm <- opop %>% filter(pid == ffmmmmm) %>% pull(mom) %>% ze_na()
  ffmmmmmf <- opop %>% filter(pid == ffmmmmm) %>% pull(pop) %>% ze_na()
  ffmmmmfm <- opop %>% filter(pid == ffmmmmf) %>% pull(mom) %>% ze_na()
  ffmmmmff <- opop %>% filter(pid == ffmmmmf) %>% pull(pop) %>% ze_na()
  ffmmmfmm <- opop %>% filter(pid == ffmmmfm) %>% pull(mom) %>% ze_na()
  ffmmmfmf <- opop %>% filter(pid == ffmmmfm) %>% pull(pop) %>% ze_na()
  ffmmmffm <- opop %>% filter(pid == ffmmmff) %>% pull(mom) %>% ze_na()
  ffmmmfff <- opop %>% filter(pid == ffmmmff) %>% pull(pop) %>% ze_na()
  ffmmfmmm <- opop %>% filter(pid == ffmmfmm) %>% pull(mom) %>% ze_na()
  ffmmfmmf <- opop %>% filter(pid == ffmmfmm) %>% pull(pop) %>% ze_na()
  ffmmfmfm <- opop %>% filter(pid == ffmmfmf) %>% pull(mom) %>% ze_na()
  ffmmfmff <- opop %>% filter(pid == ffmmfmf) %>% pull(pop) %>% ze_na()
  ffmmffmm <- opop %>% filter(pid == ffmmffm) %>% pull(mom) %>% ze_na()
  ffmmffmf <- opop %>% filter(pid == ffmmffm) %>% pull(pop) %>% ze_na()
  ffmmfffm <- opop %>% filter(pid == ffmmfff) %>% pull(mom) %>% ze_na()
  ffmmffff <- opop %>% filter(pid == ffmmfff) %>% pull(pop) %>% ze_na()
  ffmfmmmm <- opop %>% filter(pid == ffmfmmm) %>% pull(mom) %>% ze_na()
  ffmfmmmf <- opop %>% filter(pid == ffmfmmm) %>% pull(pop) %>% ze_na()
  ffmfmmfm <- opop %>% filter(pid == ffmfmmf) %>% pull(mom) %>% ze_na()
  ffmfmmff <- opop %>% filter(pid == ffmfmmf) %>% pull(pop) %>% ze_na()
  ffmfmfmm <- opop %>% filter(pid == ffmfmfm) %>% pull(mom) %>% ze_na()
  ffmfmfmf <- opop %>% filter(pid == ffmfmfm) %>% pull(pop) %>% ze_na()
  ffmfmffm <- opop %>% filter(pid == ffmfmff) %>% pull(mom) %>% ze_na()
  ffmfmfff <- opop %>% filter(pid == ffmfmff) %>% pull(pop) %>% ze_na()
  ffmffmmm <- opop %>% filter(pid == ffmffmm) %>% pull(mom) %>% ze_na()
  ffmffmmf <- opop %>% filter(pid == ffmffmm) %>% pull(pop) %>% ze_na()
  ffmffmfm <- opop %>% filter(pid == ffmffmf) %>% pull(mom) %>% ze_na()
  ffmffmff <- opop %>% filter(pid == ffmffmf) %>% pull(pop) %>% ze_na()
  ffmfffmm <- opop %>% filter(pid == ffmfffm) %>% pull(mom) %>% ze_na()
  ffmfffmf <- opop %>% filter(pid == ffmfffm) %>% pull(pop) %>% ze_na()
  ffmffffm <- opop %>% filter(pid == ffmffff) %>% pull(mom) %>% ze_na()
  ffmfffff <- opop %>% filter(pid == ffmffff) %>% pull(pop) %>% ze_na()
  fffmmmmm <- opop %>% filter(pid == fffmmmm) %>% pull(mom) %>% ze_na()
  fffmmmmf <- opop %>% filter(pid == fffmmmm) %>% pull(pop) %>% ze_na()
  fffmmmfm <- opop %>% filter(pid == fffmmmf) %>% pull(mom) %>% ze_na()
  fffmmmff <- opop %>% filter(pid == fffmmmf) %>% pull(pop) %>% ze_na()
  fffmmfmm <- opop %>% filter(pid == fffmmfm) %>% pull(mom) %>% ze_na()
  fffmmfmf <- opop %>% filter(pid == fffmmfm) %>% pull(pop) %>% ze_na()
  fffmmffm <- opop %>% filter(pid == fffmmff) %>% pull(mom) %>% ze_na()
  fffmmfff <- opop %>% filter(pid == fffmmff) %>% pull(pop) %>% ze_na()
  fffmfmmm <- opop %>% filter(pid == fffmfmm) %>% pull(mom) %>% ze_na()
  fffmfmmf <- opop %>% filter(pid == fffmfmm) %>% pull(pop) %>% ze_na()
  fffmfmfm <- opop %>% filter(pid == fffmfmf) %>% pull(mom) %>% ze_na()
  fffmfmff <- opop %>% filter(pid == fffmfmf) %>% pull(pop) %>% ze_na()
  fffmffmm <- opop %>% filter(pid == fffmffm) %>% pull(mom) %>% ze_na()
  fffmffmf <- opop %>% filter(pid == fffmffm) %>% pull(pop) %>% ze_na()
  fffmfffm <- opop %>% filter(pid == fffmfff) %>% pull(mom) %>% ze_na()
  fffmffff <- opop %>% filter(pid == fffmfff) %>% pull(pop) %>% ze_na()
  ffffmmmm <- opop %>% filter(pid == ffffmmm) %>% pull(mom) %>% ze_na()
  ffffmmmf <- opop %>% filter(pid == ffffmmm) %>% pull(pop) %>% ze_na()
  ffffmmfm <- opop %>% filter(pid == ffffmmf) %>% pull(mom) %>% ze_na()
  ffffmmff <- opop %>% filter(pid == ffffmmf) %>% pull(pop) %>% ze_na()
  ffffmfmm <- opop %>% filter(pid == ffffmfm) %>% pull(mom) %>% ze_na()
  ffffmfmf <- opop %>% filter(pid == ffffmfm) %>% pull(pop) %>% ze_na()
  ffffmffm <- opop %>% filter(pid == ffffmff) %>% pull(mom) %>% ze_na()
  ffffmfff <- opop %>% filter(pid == ffffmff) %>% pull(pop) %>% ze_na()
  fffffmmm <- opop %>% filter(pid == fffffmm) %>% pull(mom) %>% ze_na()
  fffffmmf <- opop %>% filter(pid == fffffmm) %>% pull(pop) %>% ze_na()
  fffffmfm <- opop %>% filter(pid == fffffmf) %>% pull(mom) %>% ze_na()
  fffffmff <- opop %>% filter(pid == fffffmf) %>% pull(pop) %>% ze_na()
  ffffffmm <- opop %>% filter(pid == ffffffm) %>% pull(mom) %>% ze_na()
  ffffffmf <- opop %>% filter(pid == ffffffm) %>% pull(pop) %>% ze_na()
  fffffffm <- opop %>% filter(pid == fffffff) %>% pull(mom) %>% ze_na()
  ffffffff <- opop %>% filter(pid == fffffff) %>% pull(pop) %>% ze_na()
  
  
  ## Bind all ego's direct ancestors in a tibble and include their information from .opop
  
  ancestors_opop <- tibble(kin_type = c("ego", # 1 
                                        "m", "f", # 2 
                                        "mm", "mf", "fm", "ff", # 3 
                                        "mmm", "mmf", "mfm", "mff", # 4
                                        "fmm", "fmf", "ffm", "fff", # 4
                                        "mmmm", "mmmf", "mmfm", "mmff", "mfmm", "mfmf", "mffm", "mfff", # 5
                                        "fmmm", "fmmf", "fmfm", "fmff", "ffmm", "ffmf", "fffm", "ffff", # 5
                                        "mmmmm", "mmmmf", "mmmfm", "mmmff", "mmfmm", "mmfmf", "mmffm", "mmfff", # 6
                                        "mfmmm", "mfmmf", "mfmfm", "mfmff", "mffmm", "mffmf", "mfffm", "mffff", # 6
                                        "fmmmm", "fmmmf", "fmmfm", "fmmff", "fmfmm", "fmfmf", "fmffm", "fmfff", # 6
                                        "ffmmm", "ffmmf", "ffmfm", "ffmff", "fffmm", "fffmf", "ffffm", "fffff", # 6
                                        "mmmmmm", "mmmmmf", "mmmmfm", "mmmmff", "mmmfmm", "mmmfmf", "mmmffm", "mmmfff", # 7
                                        "mmfmmm", "mmfmmf", "mmfmfm", "mmfmff", "mmffmm", "mmffmf", "mmfffm", "mmffff", # 7
                                        "mfmmmm", "mfmmmf", "mfmmfm", "mfmmff", "mfmfmm", "mfmfmf", "mfmffm", "mfmfff", # 7
                                        "mffmmm", "mffmmf", "mffmfm", "mffmff", "mfffmm", "mfffmf", "mffffm", "mfffff", # 7
                                        "fmmmmm", "fmmmmf", "fmmmfm", "fmmmff", "fmmfmm", "fmmfmf", "fmmffm", "fmmfff", # 7
                                        "fmfmmm", "fmfmmf", "fmfmfm", "fmfmff", "fmffmm", "fmffmf", "fmfffm", "fmffff", # 7
                                        "ffmmmm", "ffmmmf", "ffmmfm", "ffmmff", "ffmfmm", "ffmfmf", "ffmffm", "ffmfff", # 7
                                        "fffmmm", "fffmmf", "fffmfm", "fffmff", "ffffmm", "ffffmf", "fffffm", "ffffff",  # 7
                                        "mmmmmmm", "mmmmmmf", "mmmmmfm", "mmmmmff", "mmmmfmm", "mmmmfmf", "mmmmffm", "mmmmfff", # 8
                                        "mmmfmmm", "mmmfmmf", "mmmfmfm", "mmmfmff", "mmmffmm", "mmmffmf", "mmmfffm", "mmmffff", # 8 
                                        "mmfmmmm", "mmfmmmf", "mmfmmfm", "mmfmmff", "mmfmfmm", "mmfmfmf", "mmfmffm", "mmfmfff", # 8 
                                        "mmffmmm", "mmffmmf", "mmffmfm", "mmffmff", "mmfffmm", "mmfffmf", "mmffffm", "mmfffff", # 8 
                                        "mfmmmmm", "mfmmmmf", "mfmmmfm", "mfmmmff", "mfmmfmm", "mfmmfmf", "mfmmffm", "mfmmfff", # 8 
                                        "mfmfmmm", "mfmfmmf", "mfmfmfm", "mfmfmff", "mfmffmm", "mfmffmf", "mfmfffm", "mfmffff", # 8
                                        "mffmmmm", "mffmmmf", "mffmmfm", "mffmmff", "mffmfmm", "mffmfmf", "mffmffm", "mffmfff", # 8
                                        "mfffmmm", "mfffmmf", "mfffmfm", "mfffmff", "mffffmm", "mffffmf", "mfffffm", "mffffff", # 8
                                        "fmmmmmm", "fmmmmmf", "fmmmmfm", "fmmmmff", "fmmmfmm", "fmmmfmf", "fmmmffm", "fmmmfff", # 8
                                        "fmmfmmm", "fmmfmmf", "fmmfmfm", "fmmfmff", "fmmffmm", "fmmffmf", "fmmfffm", "fmmffff", # 8
                                        "fmfmmmm", "fmfmmmf", "fmfmmfm", "fmfmmff", "fmfmfmm", "fmfmfmf", "fmfmffm", "fmfmfff", # 8
                                        "fmffmmm", "fmffmmf", "fmffmfm", "fmffmff", "fmfffmm", "fmfffmf", "fmffffm", "fmfffff", # 8 
                                        "ffmmmmm", "ffmmmmf", "ffmmmfm", "ffmmmff", "ffmmfmm", "ffmmfmf", "ffmmffm", "ffmmfff", # 8
                                        "ffmfmmm", "ffmfmmf", "ffmfmfm", "ffmfmff", "ffmffmm", "ffmffmf", "ffmfffm", "ffmffff", # 8
                                        "fffmmmm", "fffmmmf", "fffmmfm", "fffmmff", "fffmfmm", "fffmfmf", "fffmffm", "fffmfff", # 8
                                        "ffffmmm", "ffffmmf", "ffffmfm", "ffffmff", "fffffmm", "fffffmf", "ffffffm", "fffffff",  # 8
                                        "mmmmmmmm", "mmmmmmmf", "mmmmmmfm", "mmmmmmff", "mmmmmfmm", "mmmmmfmf", "mmmmmffm", "mmmmmfff", # 9
                                        "mmmmfmmm", "mmmmfmmf", "mmmmfmfm", "mmmmfmff", "mmmmffmm", "mmmmffmf", "mmmmfffm", "mmmmffff", # 9
                                        "mmmfmmmm", "mmmfmmmf", "mmmfmmfm", "mmmfmmff", "mmmfmfmm", "mmmfmfmf", "mmmfmffm", "mmmfmfff", # 9
                                        "mmmffmmm", "mmmffmmf", "mmmffmfm", "mmmffmff", "mmmfffmm", "mmmfffmf", "mmmffffm", "mmmfffff", # 9
                                        "mmfmmmmm", "mmfmmmmf", "mmfmmmfm", "mmfmmmff", "mmfmmfmm", "mmfmmfmf", "mmfmmffm", "mmfmmfff", # 9
                                        "mmfmfmmm", "mmfmfmmf", "mmfmfmfm", "mmfmfmff", "mmfmffmm", "mmfmffmf", "mmfmfffm", "mmfmffff", # 9
                                        "mmffmmmm", "mmffmmmf", "mmffmmfm", "mmffmmff", "mmffmfmm", "mmffmfmf", "mmffmffm", "mmffmfff", # 9
                                        "mmfffmmm", "mmfffmmf", "mmfffmfm", "mmfffmff", "mmffffmm", "mmffffmf", "mmfffffm", "mmffffff", # 9
                                        "mfmmmmmm", "mfmmmmmf", "mfmmmmfm", "mfmmmmff", "mfmmmfmm", "mfmmmfmf", "mfmmmffm", "mfmmmfff", # 9
                                        "mfmmfmmm", "mfmmfmmf", "mfmmfmfm", "mfmmfmff", "mfmmffmm", "mfmmffmf", "mfmmfffm", "mfmmffff", # 9
                                        "mfmfmmmm", "mfmfmmmf", "mfmfmmfm", "mfmfmmff", "mfmfmfmm", "mfmfmfmf", "mfmfmffm", "mfmfmfff", # 9
                                        "mfmffmmm", "mfmffmmf", "mfmffmfm", "mfmffmff", "mfmfffmm", "mfmfffmf", "mfmffffm", "mfmfffff", # 9
                                        "mffmmmmm", "mffmmmmf", "mffmmmfm", "mffmmmff", "mffmmfmm", "mffmmfmf", "mffmmffm", "mffmmfff", # 9
                                        "mffmfmmm", "mffmfmmf", "mffmfmfm", "mffmfmff", "mffmffmm", "mffmffmf", "mffmfffm", "mffmffff", # 9
                                        "mfffmmmm", "mfffmmmf", "mfffmmfm", "mfffmmff", "mfffmfmm", "mfffmfmf", "mfffmffm", "mfffmfff", # 9
                                        "mffffmmm", "mffffmmf", "mffffmfm", "mffffmff", "mfffffmm", "mfffffmf", "mffffffm", "mfffffff", # 9
                                        "fmmmmmmm", "fmmmmmmf", "fmmmmmfm", "fmmmmmff", "fmmmmfmm", "fmmmmfmf", "fmmmmffm", "fmmmmfff", # 9
                                        "fmmmfmmm", "fmmmfmmf", "fmmmfmfm", "fmmmfmff", "fmmmffmm", "fmmmffmf", "fmmmfffm", "fmmmffff", # 9
                                        "fmmfmmmm", "fmmfmmmf", "fmmfmmfm", "fmmfmmff", "fmmfmfmm", "fmmfmfmf", "fmmfmffm", "fmmfmfff", # 9
                                        "fmmffmmm", "fmmffmmf", "fmmffmfm", "fmmffmff", "fmmfffmm", "fmmfffmf", "fmmffffm", "fmmfffff", # 9
                                        "fmfmmmmm", "fmfmmmmf", "fmfmmmfm", "fmfmmmff", "fmfmmfmm", "fmfmmfmf", "fmfmmffm", "fmfmmfff", # 9
                                        "fmfmfmmm", "fmfmfmmf", "fmfmfmfm", "fmfmfmff", "fmfmffmm", "fmfmffmf", "fmfmfffm", "fmfmffff", # 9
                                        "fmffmmmm", "fmffmmmf", "fmffmmfm", "fmffmmff", "fmffmfmm", "fmffmfmf", "fmffmffm", "fmffmfff", # 9
                                        "fmfffmmm", "fmfffmmf", "fmfffmfm", "fmfffmff", "fmffffmm", "fmffffmf", "fmfffffm", "fmffffff", # 9
                                        "ffmmmmmm", "ffmmmmmf", "ffmmmmfm", "ffmmmmff", "ffmmmfmm", "ffmmmfmf", "ffmmmffm", "ffmmmfff", # 9
                                        "ffmmfmmm", "ffmmfmmf", "ffmmfmfm", "ffmmfmff", "ffmmffmm", "ffmmffmf", "ffmmfffm", "ffmmffff", # 9
                                        "ffmfmmmm", "ffmfmmmf", "ffmfmmfm", "ffmfmmff", "ffmfmfmm", "ffmfmfmf", "ffmfmffm", "ffmfmfff", # 9
                                        "ffmffmmm", "ffmffmmf", "ffmffmfm", "ffmffmff", "ffmfffmm", "ffmfffmf", "ffmffffm", "ffmfffff", # 9
                                        "fffmmmmm", "fffmmmmf", "fffmmmfm", "fffmmmff", "fffmmfmm", "fffmmfmf", "fffmmffm", "fffmmfff", # 9
                                        "fffmfmmm", "fffmfmmf", "fffmfmfm", "fffmfmff", "fffmffmm", "fffmffmf", "fffmfffm", "fffmffff", # 9
                                        "ffffmmmm", "ffffmmmf", "ffffmmfm", "ffffmmff", "ffffmfmm", "ffffmfmf", "ffffmffm", "ffffmfff", # 9
                                        "fffffmmm", "fffffmmf", "fffffmfm", "fffffmff", "ffffffmm", "ffffffmf", "fffffffm", "ffffffff"  # 9
  ), 
  pid = c(ego, # 1
          m, f, # 2
          mm, mf, fm, ff, # 3
          mmm, mmf, mfm, mff, # 4
          fmm, fmf, ffm, fff, # 4
          mmmm, mmmf, mmfm, mmff, mfmm, mfmf, mffm, mfff, # 5
          fmmm, fmmf, fmfm, fmff, ffmm, ffmf, fffm, ffff, # 5
          mmmmm, mmmmf, mmmfm, mmmff, mmfmm, mmfmf, mmffm, mmfff, # 6
          mfmmm, mfmmf, mfmfm, mfmff, mffmm, mffmf, mfffm, mffff, # 6
          fmmmm, fmmmf, fmmfm, fmmff, fmfmm, fmfmf, fmffm, fmfff, # 6
          ffmmm, ffmmf, ffmfm, ffmff, fffmm, fffmf, ffffm, fffff, # 6
          mmmmmm, mmmmmf, mmmmfm, mmmmff, mmmfmm, mmmfmf, mmmffm, mmmfff, # 7
          mmfmmm, mmfmmf, mmfmfm, mmfmff, mmffmm, mmffmf, mmfffm, mmffff, # 7
          mfmmmm, mfmmmf, mfmmfm, mfmmff, mfmfmm, mfmfmf, mfmffm, mfmfff, # 7
          mffmmm, mffmmf, mffmfm, mffmff, mfffmm, mfffmf, mffffm, mfffff, # 7
          fmmmmm, fmmmmf, fmmmfm, fmmmff, fmmfmm, fmmfmf, fmmffm, fmmfff, # 7
          fmfmmm, fmfmmf, fmfmfm, fmfmff, fmffmm, fmffmf, fmfffm, fmffff, # 7
          ffmmmm, ffmmmf, ffmmfm, ffmmff, ffmfmm, ffmfmf, ffmffm, ffmfff, # 7
          fffmmm, fffmmf, fffmfm, fffmff, ffffmm, ffffmf, fffffm, ffffff,  # 7
          mmmmmmm, mmmmmmf, mmmmmfm, mmmmmff, mmmmfmm, mmmmfmf, mmmmffm, mmmmfff, # 8
          mmmfmmm, mmmfmmf, mmmfmfm, mmmfmff, mmmffmm, mmmffmf, mmmfffm, mmmffff, # 8 
          mmfmmmm, mmfmmmf, mmfmmfm, mmfmmff, mmfmfmm, mmfmfmf, mmfmffm, mmfmfff, # 8 
          mmffmmm, mmffmmf, mmffmfm, mmffmff, mmfffmm, mmfffmf, mmffffm, mmfffff, # 8 
          mfmmmmm, mfmmmmf, mfmmmfm, mfmmmff, mfmmfmm, mfmmfmf, mfmmffm, mfmmfff, # 8 
          mfmfmmm, mfmfmmf, mfmfmfm, mfmfmff, mfmffmm, mfmffmf, mfmfffm, mfmffff, # 8
          mffmmmm, mffmmmf, mffmmfm, mffmmff, mffmfmm, mffmfmf, mffmffm, mffmfff, # 8
          mfffmmm, mfffmmf, mfffmfm, mfffmff, mffffmm, mffffmf, mfffffm, mffffff, # 8
          fmmmmmm, fmmmmmf, fmmmmfm, fmmmmff, fmmmfmm, fmmmfmf, fmmmffm, fmmmfff, # 8
          fmmfmmm, fmmfmmf, fmmfmfm, fmmfmff, fmmffmm, fmmffmf, fmmfffm, fmmffff, # 8
          fmfmmmm, fmfmmmf, fmfmmfm, fmfmmff, fmfmfmm, fmfmfmf, fmfmffm, fmfmfff, # 8
          fmffmmm, fmffmmf, fmffmfm, fmffmff, fmfffmm, fmfffmf, fmffffm, fmfffff, # 8 
          ffmmmmm, ffmmmmf, ffmmmfm, ffmmmff, ffmmfmm, ffmmfmf, ffmmffm, ffmmfff, # 8
          ffmfmmm, ffmfmmf, ffmfmfm, ffmfmff, ffmffmm, ffmffmf, ffmfffm, ffmffff, # 8
          fffmmmm, fffmmmf, fffmmfm, fffmmff, fffmfmm, fffmfmf, fffmffm, fffmfff, # 8
          ffffmmm, ffffmmf, ffffmfm, ffffmff, fffffmm, fffffmf, ffffffm, fffffff,  # 8
          mmmmmmmm, mmmmmmmf, mmmmmmfm, mmmmmmff, mmmmmfmm, mmmmmfmf, mmmmmffm, mmmmmfff, # 9
          mmmmfmmm, mmmmfmmf, mmmmfmfm, mmmmfmff, mmmmffmm, mmmmffmf, mmmmfffm, mmmmffff, # 9
          mmmfmmmm, mmmfmmmf, mmmfmmfm, mmmfmmff, mmmfmfmm, mmmfmfmf, mmmfmffm, mmmfmfff, # 9
          mmmffmmm, mmmffmmf, mmmffmfm, mmmffmff, mmmfffmm, mmmfffmf, mmmffffm, mmmfffff, # 9
          mmfmmmmm, mmfmmmmf, mmfmmmfm, mmfmmmff, mmfmmfmm, mmfmmfmf, mmfmmffm, mmfmmfff, # 9
          mmfmfmmm, mmfmfmmf, mmfmfmfm, mmfmfmff, mmfmffmm, mmfmffmf, mmfmfffm, mmfmffff, # 9
          mmffmmmm, mmffmmmf, mmffmmfm, mmffmmff, mmffmfmm, mmffmfmf, mmffmffm, mmffmfff, # 9
          mmfffmmm, mmfffmmf, mmfffmfm, mmfffmff, mmffffmm, mmffffmf, mmfffffm, mmffffff, # 9
          mfmmmmmm, mfmmmmmf, mfmmmmfm, mfmmmmff, mfmmmfmm, mfmmmfmf, mfmmmffm, mfmmmfff, # 9
          mfmmfmmm, mfmmfmmf, mfmmfmfm, mfmmfmff, mfmmffmm, mfmmffmf, mfmmfffm, mfmmffff, # 9
          mfmfmmmm, mfmfmmmf, mfmfmmfm, mfmfmmff, mfmfmfmm, mfmfmfmf, mfmfmffm, mfmfmfff, # 9
          mfmffmmm, mfmffmmf, mfmffmfm, mfmffmff, mfmfffmm, mfmfffmf, mfmffffm, mfmfffff, # 9
          mffmmmmm, mffmmmmf, mffmmmfm, mffmmmff, mffmmfmm, mffmmfmf, mffmmffm, mffmmfff, # 9
          mffmfmmm, mffmfmmf, mffmfmfm, mffmfmff, mffmffmm, mffmffmf, mffmfffm, mffmffff, # 9
          mfffmmmm, mfffmmmf, mfffmmfm, mfffmmff, mfffmfmm, mfffmfmf, mfffmffm, mfffmfff, # 9
          mffffmmm, mffffmmf, mffffmfm, mffffmff, mfffffmm, mfffffmf, mffffffm, mfffffff, # 9
          fmmmmmmm, fmmmmmmf, fmmmmmfm, fmmmmmff, fmmmmfmm, fmmmmfmf, fmmmmffm, fmmmmfff, # 9
          fmmmfmmm, fmmmfmmf, fmmmfmfm, fmmmfmff, fmmmffmm, fmmmffmf, fmmmfffm, fmmmffff, # 9
          fmmfmmmm, fmmfmmmf, fmmfmmfm, fmmfmmff, fmmfmfmm, fmmfmfmf, fmmfmffm, fmmfmfff, # 9
          fmmffmmm, fmmffmmf, fmmffmfm, fmmffmff, fmmfffmm, fmmfffmf, fmmffffm, fmmfffff, # 9
          fmfmmmmm, fmfmmmmf, fmfmmmfm, fmfmmmff, fmfmmfmm, fmfmmfmf, fmfmmffm, fmfmmfff, # 9
          fmfmfmmm, fmfmfmmf, fmfmfmfm, fmfmfmff, fmfmffmm, fmfmffmf, fmfmfffm, fmfmffff, # 9
          fmffmmmm, fmffmmmf, fmffmmfm, fmffmmff, fmffmfmm, fmffmfmf, fmffmffm, fmffmfff, # 9
          fmfffmmm, fmfffmmf, fmfffmfm, fmfffmff, fmffffmm, fmffffmf, fmfffffm, fmffffff, # 9
          ffmmmmmm, ffmmmmmf, ffmmmmfm, ffmmmmff, ffmmmfmm, ffmmmfmf, ffmmmffm, ffmmmfff, # 9
          ffmmfmmm, ffmmfmmf, ffmmfmfm, ffmmfmff, ffmmffmm, ffmmffmf, ffmmfffm, ffmmffff, # 9
          ffmfmmmm, ffmfmmmf, ffmfmmfm, ffmfmmff, ffmfmfmm, ffmfmfmf, ffmfmffm, ffmfmfff, # 9
          ffmffmmm, ffmffmmf, ffmffmfm, ffmffmff, ffmfffmm, ffmfffmf, ffmffffm, ffmfffff, # 9
          fffmmmmm, fffmmmmf, fffmmmfm, fffmmmff, fffmmfmm, fffmmfmf, fffmmffm, fffmmfff, # 9
          fffmfmmm, fffmfmmf, fffmfmfm, fffmfmff, fffmffmm, fffmffmf, fffmfffm, fffmffff, # 9
          ffffmmmm, ffffmmmf, ffffmmfm, ffffmmff, ffffmfmm, ffffmfmf, ffffmffm, ffffmfff, # 9
          fffffmmm, fffffmmf, fffffmfm, fffffmff, ffffffmm, ffffffmf, fffffffm, ffffffff  # 9
  ),
  ego_id = rep(ego, length(pid))) %>% 
    filter(!is.na(pid)) %>% 
    distinct(pid, .keep_all = TRUE) %>%
    left_join(select(opop, c(pid, fem, dob, dod, mom, marid, mstat)),
              by = "pid")
  
  return(ancestors_opop)
  
}