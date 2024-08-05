#----------------------------------------------------------------------------------------------------
# Function to retrieve pids of direct ancestors up to the 9th generation of a given ego(s)
# from SOCSIM microsimulation outputs

# Created on 23-09-2022
# Last modified on 14-06-2023

# This is a modified version of the script "U:/SOCSIM/SOCSIM_Genealogies/SWE_PAA/R_old/Function_Ancestors_Old.R" 
#----------------------------------------------------------------------------------------------------

retrieve_ancestors <- function(egos, opop = opop) {
  
  opop2 <- opop %>% 
    select(pid, mom, pop) %>% 
    # Change 0s (nothing) for the pids of the initial population
    mutate(mom = ifelse(mom == 0, NA, mom), 
           pop = ifelse(pop == 0, NA, pop))
  
  egos <- egos

  ### 3. Grandmothers
  opop2$mm <- opop2$mom[match(opop2$mom, opop2$pid)]
  opop2$mf <- opop2$mom[match(opop2$pop, opop2$pid)]
  ### 3. Grandfathers
  opop2$fm <- opop2$pop[match(opop2$mom, opop2$pid)]  
  opop2$ff <- opop2$pop[match(opop2$pop, opop2$pid)]

  ### 4. Great-grandmothers
  opop2$mmm <- opop2$mm[match(opop2$mom, opop2$pid)]
  opop2$mfm <- opop2$mf[match(opop2$mom, opop2$pid)]
  opop2$mmf <- opop2$mm[match(opop2$pop, opop2$pid)]
  opop2$mff <- opop2$mf[match(opop2$pop, opop2$pid)]
  ### 4. Great-grandfathers
  opop2$fmm <- opop2$fm[match(opop2$mom, opop2$pid)]
  opop2$ffm <- opop2$ff[match(opop2$mom, opop2$pid)] 
  opop2$fmf <- opop2$fm[match(opop2$pop, opop2$pid)]
  opop2$fff <- opop2$ff[match(opop2$pop, opop2$pid)]
  
  
  ### 5. Great-great-grandmothers
  opop2$mmmm <- opop2$mmm[match(opop2$mom, opop2$pid)]
  opop2$mmfm <- opop2$mmf[match(opop2$mom, opop2$pid)]
  opop2$mfmm <- opop2$mfm[match(opop2$mom, opop2$pid)]
  opop2$mffm <- opop2$mff[match(opop2$mom, opop2$pid)]
  opop2$mmmf <- opop2$mmm[match(opop2$pop, opop2$pid)]
  opop2$mmff <- opop2$mmf[match(opop2$pop, opop2$pid)]
  opop2$mfmf <- opop2$mfm[match(opop2$pop, opop2$pid)]
  opop2$mfff <- opop2$mff[match(opop2$pop, opop2$pid)]
  ### 5. Great-great-grandfathers
  opop2$fmmm <- opop2$fmm[match(opop2$mom, opop2$pid)]
  opop2$fmfm <- opop2$fmf[match(opop2$mom, opop2$pid)]
  opop2$ffmm <- opop2$ffm[match(opop2$mom, opop2$pid)]
  opop2$fffm <- opop2$fff[match(opop2$mom, opop2$pid)]
  opop2$fmmf <- opop2$fmm[match(opop2$pop, opop2$pid)]
  opop2$fmff <- opop2$fmf[match(opop2$pop, opop2$pid)]
  opop2$ffmf <- opop2$ffm[match(opop2$pop, opop2$pid)]
  opop2$ffff <- opop2$fff[match(opop2$pop, opop2$pid)]

  ### 6. Great-great-great-grandmothers
  opop2$mmmmm <- opop2$mmmm[match(opop2$mom, opop2$pid)]
  opop2$mmmfm <- opop2$mmmf[match(opop2$mom, opop2$pid)]
  opop2$mmfmm <- opop2$mmfm[match(opop2$mom, opop2$pid)]
  opop2$mmffm <- opop2$mmff[match(opop2$mom, opop2$pid)]
  opop2$mfmmm <- opop2$mfmm[match(opop2$mom, opop2$pid)]
  opop2$mfmfm <- opop2$mfmf[match(opop2$mom, opop2$pid)]
  opop2$mffmm <- opop2$mffm[match(opop2$mom, opop2$pid)]
  opop2$mfffm <- opop2$mfff[match(opop2$mom, opop2$pid)]
  opop2$mmmmf <- opop2$mmmm[match(opop2$pop, opop2$pid)]
  opop2$mmmff <- opop2$mmmf[match(opop2$pop, opop2$pid)]
  opop2$mmfmf <- opop2$mmfm[match(opop2$pop, opop2$pid)]
  opop2$mmfff <- opop2$mmff[match(opop2$pop, opop2$pid)]
  opop2$mfmmf <- opop2$mfmm[match(opop2$pop, opop2$pid)]
  opop2$mfmff <- opop2$mfmf[match(opop2$pop, opop2$pid)]
  opop2$mffmf <- opop2$mffm[match(opop2$pop, opop2$pid)]
  opop2$mffff <- opop2$mfff[match(opop2$pop, opop2$pid)]
  ### 6. Great-great-great-grandfathers
  opop2$fmmmm <- opop2$fmmm[match(opop2$mom, opop2$pid)]
  opop2$fmmfm <- opop2$fmmf[match(opop2$mom, opop2$pid)]
  opop2$fmfmm <- opop2$fmfm[match(opop2$mom, opop2$pid)]
  opop2$fmffm <- opop2$fmff[match(opop2$mom, opop2$pid)]
  opop2$ffmmm <- opop2$ffmm[match(opop2$mom, opop2$pid)]
  opop2$ffmfm <- opop2$ffmf[match(opop2$mom, opop2$pid)]
  opop2$fffmm <- opop2$fffm[match(opop2$mom, opop2$pid)]
  opop2$ffffm <- opop2$ffff[match(opop2$mom, opop2$pid)]
  opop2$fmmmf <- opop2$fmmm[match(opop2$pop, opop2$pid)]
  opop2$fmmff <- opop2$fmmf[match(opop2$pop, opop2$pid)]
  opop2$fmfmf <- opop2$fmfm[match(opop2$pop, opop2$pid)]
  opop2$fmfff <- opop2$fmff[match(opop2$pop, opop2$pid)]
  opop2$ffmmf <- opop2$ffmm[match(opop2$pop, opop2$pid)]
  opop2$ffmff <- opop2$ffmf[match(opop2$pop, opop2$pid)]
  opop2$fffmf <- opop2$fffm[match(opop2$pop, opop2$pid)]
  opop2$fffff <- opop2$ffff[match(opop2$pop, opop2$pid)]
  
  ### 7. Great-great-great-great-grandmothers
  opop2$mmmmmm <- opop2$mmmmm[match(opop2$mom, opop2$pid)]
  opop2$mmmmfm <- opop2$mmmmf[match(opop2$mom, opop2$pid)]
  opop2$mmmfmm <- opop2$mmmfm[match(opop2$mom, opop2$pid)]
  opop2$mmmffm <- opop2$mmmff[match(opop2$mom, opop2$pid)]
  opop2$mmfmmm <- opop2$mmfmm[match(opop2$mom, opop2$pid)]
  opop2$mmfmfm <- opop2$mmfmf[match(opop2$mom, opop2$pid)]
  opop2$mmffmm <- opop2$mmffm[match(opop2$mom, opop2$pid)]
  opop2$mmfffm <- opop2$mmfff[match(opop2$mom, opop2$pid)]
  opop2$mfmmmm <- opop2$mfmmm[match(opop2$mom, opop2$pid)]
  opop2$mfmmfm <- opop2$mfmmf[match(opop2$mom, opop2$pid)]
  opop2$mfmfmm <- opop2$mfmfm[match(opop2$mom, opop2$pid)]
  opop2$mfmffm <- opop2$mfmff[match(opop2$mom, opop2$pid)]
  opop2$mffmmm <- opop2$mffmm[match(opop2$mom, opop2$pid)]
  opop2$mffmfm <- opop2$mffmf[match(opop2$mom, opop2$pid)]
  opop2$mfffmm <- opop2$mfffm[match(opop2$mom, opop2$pid)]
  opop2$mffffm <- opop2$mffff[match(opop2$mom, opop2$pid)]
  opop2$mmmmmf <- opop2$mmmmm[match(opop2$pop, opop2$pid)]
  opop2$mmmmff <- opop2$mmmmf[match(opop2$pop, opop2$pid)]
  opop2$mmmfmf <- opop2$mmmfm[match(opop2$pop, opop2$pid)]
  opop2$mmmfff <- opop2$mmmff[match(opop2$pop, opop2$pid)]
  opop2$mmfmmf <- opop2$mmfmm[match(opop2$pop, opop2$pid)]
  opop2$mmfmff <- opop2$mmfmf[match(opop2$pop, opop2$pid)]
  opop2$mmffmf <- opop2$mmffm[match(opop2$pop, opop2$pid)]
  opop2$mmffff <- opop2$mmfff[match(opop2$pop, opop2$pid)]
  opop2$mfmmmf <- opop2$mfmmm[match(opop2$pop, opop2$pid)]
  opop2$mfmmff <- opop2$mfmmf[match(opop2$pop, opop2$pid)]
  opop2$mfmfmf <- opop2$mfmfm[match(opop2$pop, opop2$pid)]
  opop2$mfmfff <- opop2$mfmff[match(opop2$pop, opop2$pid)]
  opop2$mffmmf <- opop2$mffmm[match(opop2$pop, opop2$pid)]
  opop2$mffmff <- opop2$mffmf[match(opop2$pop, opop2$pid)]
  opop2$mfffmf <- opop2$mfffm[match(opop2$pop, opop2$pid)]
  opop2$mfffff <- opop2$mffff[match(opop2$pop, opop2$pid)]
  ### 7. Great-great-great-great-grandfathers
  opop2$fmmmmm <- opop2$fmmmm[match(opop2$mom, opop2$pid)]
  opop2$fmmmfm <- opop2$fmmmf[match(opop2$mom, opop2$pid)]
  opop2$fmmfmm <- opop2$fmmfm[match(opop2$mom, opop2$pid)]
  opop2$fmmffm <- opop2$fmmff[match(opop2$mom, opop2$pid)]
  opop2$fmfmmm <- opop2$fmfmm[match(opop2$mom, opop2$pid)]
  opop2$fmfmfm <- opop2$fmfmf[match(opop2$mom, opop2$pid)]
  opop2$fmffmm <- opop2$fmffm[match(opop2$mom, opop2$pid)]
  opop2$fmfffm <- opop2$fmfff[match(opop2$mom, opop2$pid)]
  opop2$ffmmmm <- opop2$ffmmm[match(opop2$mom, opop2$pid)]
  opop2$ffmmfm <- opop2$ffmmf[match(opop2$mom, opop2$pid)]
  opop2$ffmfmm <- opop2$ffmfm[match(opop2$mom, opop2$pid)]
  opop2$ffmffm <- opop2$ffmff[match(opop2$mom, opop2$pid)]
  opop2$fffmmm <- opop2$fffmm[match(opop2$mom, opop2$pid)]
  opop2$fffmfm <- opop2$fffmf[match(opop2$mom, opop2$pid)]
  opop2$ffffmm <- opop2$ffffm[match(opop2$mom, opop2$pid)]
  opop2$fffffm <- opop2$fffff[match(opop2$mom, opop2$pid)]
  opop2$fmmmmf <- opop2$fmmmm[match(opop2$pop, opop2$pid)]
  opop2$fmmmff <- opop2$fmmmf[match(opop2$pop, opop2$pid)]
  opop2$fmmfmf <- opop2$fmmfm[match(opop2$pop, opop2$pid)]
  opop2$fmmfff <- opop2$fmmff[match(opop2$pop, opop2$pid)]
  opop2$fmfmmf <- opop2$fmfmm[match(opop2$pop, opop2$pid)]
  opop2$fmfmff <- opop2$fmfmf[match(opop2$pop, opop2$pid)]
  opop2$fmffmf <- opop2$fmffm[match(opop2$pop, opop2$pid)]
  opop2$fmffff <- opop2$fmfff[match(opop2$pop, opop2$pid)]
  opop2$ffmmmf <- opop2$ffmmm[match(opop2$pop, opop2$pid)]
  opop2$ffmmff <- opop2$ffmmf[match(opop2$pop, opop2$pid)]
  opop2$ffmfmf <- opop2$ffmfm[match(opop2$pop, opop2$pid)]
  opop2$ffmfff <- opop2$ffmff[match(opop2$pop, opop2$pid)]
  opop2$fffmmf <- opop2$fffmm[match(opop2$pop, opop2$pid)]
  opop2$fffmff <- opop2$fffmf[match(opop2$pop, opop2$pid)]
  opop2$ffffmf <- opop2$ffffm[match(opop2$pop, opop2$pid)]
  opop2$ffffff <- opop2$fffff[match(opop2$pop, opop2$pid)]
  
  
  ### 8. Great-great-great-great-great-grandmothers
  opop2$mmmmmmm <- opop2$mmmmmm[match(opop2$mom, opop2$pid)]
  opop2$mmmmmfm <- opop2$mmmmmf[match(opop2$mom, opop2$pid)]
  opop2$mmmmfmm <- opop2$mmmmfm[match(opop2$mom, opop2$pid)]
  opop2$mmmmffm <- opop2$mmmmff[match(opop2$mom, opop2$pid)]
  opop2$mmmfmmm <- opop2$mmmfmm[match(opop2$mom, opop2$pid)]
  opop2$mmmfmfm <- opop2$mmmfmf[match(opop2$mom, opop2$pid)]
  opop2$mmmffmm <- opop2$mmmffm[match(opop2$mom, opop2$pid)]
  opop2$mmmfffm <- opop2$mmmfff[match(opop2$mom, opop2$pid)]
  opop2$mmfmmmm <- opop2$mmfmmm[match(opop2$mom, opop2$pid)]
  opop2$mmfmmfm <- opop2$mmfmmf[match(opop2$mom, opop2$pid)]
  opop2$mmfmfmm <- opop2$mmfmfm[match(opop2$mom, opop2$pid)]
  opop2$mmfmffm <- opop2$mmfmff[match(opop2$mom, opop2$pid)]
  opop2$mmffmmm <- opop2$mmffmm[match(opop2$mom, opop2$pid)]
  opop2$mmffmfm <- opop2$mmffmf[match(opop2$mom, opop2$pid)]
  opop2$mmfffmm <- opop2$mmfffm[match(opop2$mom, opop2$pid)]
  opop2$mmffffm <- opop2$mmffff[match(opop2$mom, opop2$pid)]
  opop2$mfmmmmm <- opop2$mfmmmm[match(opop2$mom, opop2$pid)]
  opop2$mfmmmfm <- opop2$mfmmmf[match(opop2$mom, opop2$pid)]
  opop2$mfmmfmm <- opop2$mfmmfm[match(opop2$mom, opop2$pid)]
  opop2$mfmmffm <- opop2$mfmmff[match(opop2$mom, opop2$pid)]
  opop2$mfmfmmm <- opop2$mfmfmm[match(opop2$mom, opop2$pid)]
  opop2$mfmfmfm <- opop2$mfmfmf[match(opop2$mom, opop2$pid)]
  opop2$mfmffmm <- opop2$mfmffm[match(opop2$mom, opop2$pid)]
  opop2$mfmfffm <- opop2$mfmfff[match(opop2$mom, opop2$pid)]
  opop2$mffmmmm <- opop2$mffmmm[match(opop2$mom, opop2$pid)]
  opop2$mffmmfm <- opop2$mffmmf[match(opop2$mom, opop2$pid)]
  opop2$mffmfmm <- opop2$mffmfm[match(opop2$mom, opop2$pid)]
  opop2$mffmffm <- opop2$mffmff[match(opop2$mom, opop2$pid)]
  opop2$mfffmmm <- opop2$mfffmm[match(opop2$mom, opop2$pid)]
  opop2$mfffmfm <- opop2$mfffmf[match(opop2$mom, opop2$pid)]
  opop2$mffffmm <- opop2$mffffm[match(opop2$mom, opop2$pid)]
  opop2$mfffffm <- opop2$mfffff[match(opop2$mom, opop2$pid)]
  opop2$mmmmmmf <- opop2$mmmmmm[match(opop2$pop, opop2$pid)]
  opop2$mmmmmff <- opop2$mmmmmf[match(opop2$pop, opop2$pid)]
  opop2$mmmmfmf <- opop2$mmmmfm[match(opop2$pop, opop2$pid)]
  opop2$mmmmfff <- opop2$mmmmff[match(opop2$pop, opop2$pid)]
  opop2$mmmfmmf <- opop2$mmmfmm[match(opop2$pop, opop2$pid)]
  opop2$mmmfmff <- opop2$mmmfmf[match(opop2$pop, opop2$pid)]
  opop2$mmmffmf <- opop2$mmmffm[match(opop2$pop, opop2$pid)]
  opop2$mmmffff <- opop2$mmmfff[match(opop2$pop, opop2$pid)]
  opop2$mmfmmmf <- opop2$mmfmmm[match(opop2$pop, opop2$pid)]
  opop2$mmfmmff <- opop2$mmfmmf[match(opop2$pop, opop2$pid)]
  opop2$mmfmfmf <- opop2$mmfmfm[match(opop2$pop, opop2$pid)]
  opop2$mmfmfff <- opop2$mmfmff[match(opop2$pop, opop2$pid)]
  opop2$mmffmmf <- opop2$mmffmm[match(opop2$pop, opop2$pid)]
  opop2$mmffmff <- opop2$mmffmf[match(opop2$pop, opop2$pid)]
  opop2$mmfffmf <- opop2$mmfffm[match(opop2$pop, opop2$pid)]
  opop2$mmfffff <- opop2$mmffff[match(opop2$pop, opop2$pid)]
  opop2$mfmmmmf <- opop2$mfmmmm[match(opop2$pop, opop2$pid)]
  opop2$mfmmmff <- opop2$mfmmmf[match(opop2$pop, opop2$pid)]
  opop2$mfmmfmf <- opop2$mfmmfm[match(opop2$pop, opop2$pid)]
  opop2$mfmmfff <- opop2$mfmmff[match(opop2$pop, opop2$pid)]
  opop2$mfmfmmf <- opop2$mfmfmm[match(opop2$pop, opop2$pid)]
  opop2$mfmfmff <- opop2$mfmfmf[match(opop2$pop, opop2$pid)]
  opop2$mfmffmf <- opop2$mfmffm[match(opop2$pop, opop2$pid)]
  opop2$mfmffff <- opop2$mfmfff[match(opop2$pop, opop2$pid)]
  opop2$mffmmmf <- opop2$mffmmm[match(opop2$pop, opop2$pid)]
  opop2$mffmmff <- opop2$mffmmf[match(opop2$pop, opop2$pid)]
  opop2$mffmfmf <- opop2$mffmfm[match(opop2$pop, opop2$pid)]
  opop2$mffmfff <- opop2$mffmff[match(opop2$pop, opop2$pid)]
  opop2$mfffmmf <- opop2$mfffmm[match(opop2$pop, opop2$pid)]
  opop2$mfffmff <- opop2$mfffmf[match(opop2$pop, opop2$pid)]
  opop2$mffffmf <- opop2$mffffm[match(opop2$pop, opop2$pid)]
  opop2$mffffff <- opop2$mfffff[match(opop2$pop, opop2$pid)]
 
  ### 8. Great-great-great-great-great-grandfathers
  opop2$fmmmmmm <- opop2$fmmmmm[match(opop2$mom, opop2$pid)]
  opop2$fmmmmfm <- opop2$fmmmmf[match(opop2$mom, opop2$pid)]
  opop2$fmmmfmm <- opop2$fmmmfm[match(opop2$mom, opop2$pid)]
  opop2$fmmmffm <- opop2$fmmmff[match(opop2$mom, opop2$pid)]
  opop2$fmmfmmm <- opop2$fmmfmm[match(opop2$mom, opop2$pid)]
  opop2$fmmfmfm <- opop2$fmmfmf[match(opop2$mom, opop2$pid)]
  opop2$fmmffmm <- opop2$fmmffm[match(opop2$mom, opop2$pid)]
  opop2$fmmfffm <- opop2$fmmfff[match(opop2$mom, opop2$pid)]
  opop2$fmfmmmm <- opop2$fmfmmm[match(opop2$mom, opop2$pid)]
  opop2$fmfmmfm <- opop2$fmfmmf[match(opop2$mom, opop2$pid)]
  opop2$fmfmfmm <- opop2$fmfmfm[match(opop2$mom, opop2$pid)]
  opop2$fmfmffm <- opop2$fmfmff[match(opop2$mom, opop2$pid)]
  opop2$fmffmmm <- opop2$fmffmm[match(opop2$mom, opop2$pid)]
  opop2$fmffmfm <- opop2$fmffmf[match(opop2$mom, opop2$pid)]
  opop2$fmfffmm <- opop2$fmfffm[match(opop2$mom, opop2$pid)]
  opop2$fmffffm <- opop2$fmffff[match(opop2$mom, opop2$pid)]
  opop2$ffmmmmm <- opop2$ffmmmm[match(opop2$mom, opop2$pid)]
  opop2$ffmmmfm <- opop2$ffmmmf[match(opop2$mom, opop2$pid)]
  opop2$ffmmfmm <- opop2$ffmmfm[match(opop2$mom, opop2$pid)]
  opop2$ffmmffm <- opop2$ffmmff[match(opop2$mom, opop2$pid)]
  opop2$ffmfmmm <- opop2$ffmfmm[match(opop2$mom, opop2$pid)]
  opop2$ffmfmfm <- opop2$ffmfmf[match(opop2$mom, opop2$pid)]
  opop2$ffmffmm <- opop2$ffmffm[match(opop2$mom, opop2$pid)]
  opop2$ffmfffm <- opop2$ffmfff[match(opop2$mom, opop2$pid)]
  opop2$fffmmmm <- opop2$fffmmm[match(opop2$mom, opop2$pid)]
  opop2$fffmmfm <- opop2$fffmmf[match(opop2$mom, opop2$pid)]
  opop2$fffmfmm <- opop2$fffmfm[match(opop2$mom, opop2$pid)]
  opop2$fffmffm <- opop2$fffmff[match(opop2$mom, opop2$pid)]
  opop2$ffffmmm <- opop2$ffffmm[match(opop2$mom, opop2$pid)]
  opop2$ffffmfm <- opop2$ffffmf[match(opop2$mom, opop2$pid)]
  opop2$fffffmm <- opop2$fffffm[match(opop2$mom, opop2$pid)]
  opop2$ffffffm <- opop2$ffffff[match(opop2$mom, opop2$pid)]
  opop2$fmmmmmf <- opop2$fmmmmm[match(opop2$pop, opop2$pid)]
  opop2$fmmmmff <- opop2$fmmmmf[match(opop2$pop, opop2$pid)]
  opop2$fmmmfmf <- opop2$fmmmfm[match(opop2$pop, opop2$pid)]
  opop2$fmmmfff <- opop2$fmmmff[match(opop2$pop, opop2$pid)]
  opop2$fmmfmmf <- opop2$fmmfmm[match(opop2$pop, opop2$pid)]
  opop2$fmmfmff <- opop2$fmmfmf[match(opop2$pop, opop2$pid)]
  opop2$fmmffmf <- opop2$fmmffm[match(opop2$pop, opop2$pid)]
  opop2$fmmffff <- opop2$fmmfff[match(opop2$pop, opop2$pid)]
  opop2$fmfmmmf <- opop2$fmfmmm[match(opop2$pop, opop2$pid)]
  opop2$fmfmmff <- opop2$fmfmmf[match(opop2$pop, opop2$pid)]
  opop2$fmfmfmf <- opop2$fmfmfm[match(opop2$pop, opop2$pid)]
  opop2$fmfmfff <- opop2$fmfmff[match(opop2$pop, opop2$pid)]
  opop2$fmffmmf <- opop2$fmffmm[match(opop2$pop, opop2$pid)]
  opop2$fmffmff <- opop2$fmffmf[match(opop2$pop, opop2$pid)]
  opop2$fmfffmf <- opop2$fmfffm[match(opop2$pop, opop2$pid)]
  opop2$fmfffff <- opop2$fmffff[match(opop2$pop, opop2$pid)]
  opop2$ffmmmmf <- opop2$ffmmmm[match(opop2$pop, opop2$pid)]
  opop2$ffmmmff <- opop2$ffmmmf[match(opop2$pop, opop2$pid)]
  opop2$ffmmfmf <- opop2$ffmmfm[match(opop2$pop, opop2$pid)]
  opop2$ffmmfff <- opop2$ffmmff[match(opop2$pop, opop2$pid)]
  opop2$ffmfmmf <- opop2$ffmfmm[match(opop2$pop, opop2$pid)]
  opop2$ffmfmff <- opop2$ffmfmf[match(opop2$pop, opop2$pid)]
  opop2$ffmffmf <- opop2$ffmffm[match(opop2$pop, opop2$pid)]
  opop2$ffmffff <- opop2$ffmfff[match(opop2$pop, opop2$pid)]
  opop2$fffmmmf <- opop2$fffmmm[match(opop2$pop, opop2$pid)]
  opop2$fffmmff <- opop2$fffmmf[match(opop2$pop, opop2$pid)]
  opop2$fffmfmf <- opop2$fffmfm[match(opop2$pop, opop2$pid)]
  opop2$fffmfff <- opop2$fffmff[match(opop2$pop, opop2$pid)]
  opop2$ffffmmf <- opop2$ffffmm[match(opop2$pop, opop2$pid)]
  opop2$ffffmff <- opop2$ffffmf[match(opop2$pop, opop2$pid)]
  opop2$fffffmf <- opop2$fffffm[match(opop2$pop, opop2$pid)]
  opop2$fffffff <- opop2$ffffff[match(opop2$pop, opop2$pid)]
  
  
  ### 9. Great-great-great-great-great-great-grandmothers
  opop2$mmmmmmmm <- opop2$mmmmmmm[match(opop2$mom, opop2$pid)]
  opop2$mmmmmmfm <- opop2$mmmmmmf[match(opop2$mom, opop2$pid)]
  opop2$mmmmmfmm <- opop2$mmmmmfm[match(opop2$mom, opop2$pid)]
  opop2$mmmmmffm <- opop2$mmmmmff[match(opop2$mom, opop2$pid)]
  opop2$mmmmfmmm <- opop2$mmmmfmm[match(opop2$mom, opop2$pid)]
  opop2$mmmmfmfm <- opop2$mmmmfmf[match(opop2$mom, opop2$pid)]
  opop2$mmmmffmm <- opop2$mmmmffm[match(opop2$mom, opop2$pid)]
  opop2$mmmmfffm <- opop2$mmmmfff[match(opop2$mom, opop2$pid)]
  opop2$mmmfmmmm <- opop2$mmmfmmm[match(opop2$mom, opop2$pid)]
  opop2$mmmfmmfm <- opop2$mmmfmmf[match(opop2$mom, opop2$pid)]
  opop2$mmmfmfmm <- opop2$mmmfmfm[match(opop2$mom, opop2$pid)]
  opop2$mmmfmffm <- opop2$mmmfmff[match(opop2$mom, opop2$pid)]
  opop2$mmmffmmm <- opop2$mmmffmm[match(opop2$mom, opop2$pid)]
  opop2$mmmffmfm <- opop2$mmmffmf[match(opop2$mom, opop2$pid)]
  opop2$mmmfffmm <- opop2$mmmfffm[match(opop2$mom, opop2$pid)]
  opop2$mmmffffm <- opop2$mmmffff[match(opop2$mom, opop2$pid)]
  opop2$mmfmmmmm <- opop2$mmfmmmm[match(opop2$mom, opop2$pid)]
  opop2$mmfmmmfm <- opop2$mmfmmmf[match(opop2$mom, opop2$pid)]
  opop2$mmfmmfmm <- opop2$mmfmmfm[match(opop2$mom, opop2$pid)]
  opop2$mmfmmffm <- opop2$mmfmmff[match(opop2$mom, opop2$pid)]
  opop2$mmfmfmmm <- opop2$mmfmfmm[match(opop2$mom, opop2$pid)]
  opop2$mmfmfmfm <- opop2$mmfmfmf[match(opop2$mom, opop2$pid)]
  opop2$mmfmffmm <- opop2$mmfmffm[match(opop2$mom, opop2$pid)]
  opop2$mmfmfffm <- opop2$mmfmfff[match(opop2$mom, opop2$pid)]
  opop2$mmffmmmm <- opop2$mmffmmm[match(opop2$mom, opop2$pid)]
  opop2$mmffmmfm <- opop2$mmffmmf[match(opop2$mom, opop2$pid)]
  opop2$mmffmfmm <- opop2$mmffmfm[match(opop2$mom, opop2$pid)]
  opop2$mmffmffm <- opop2$mmffmff[match(opop2$mom, opop2$pid)]
  opop2$mmfffmmm <- opop2$mmfffmm[match(opop2$mom, opop2$pid)]
  opop2$mmfffmfm <- opop2$mmfffmf[match(opop2$mom, opop2$pid)]
  opop2$mmffffmm <- opop2$mmffffm[match(opop2$mom, opop2$pid)]
  opop2$mmfffffm <- opop2$mmfffff[match(opop2$mom, opop2$pid)]
  opop2$mfmmmmmm <- opop2$mfmmmmm[match(opop2$mom, opop2$pid)]
  opop2$mfmmmmfm <- opop2$mfmmmmf[match(opop2$mom, opop2$pid)]
  opop2$mfmmmfmm <- opop2$mfmmmfm[match(opop2$mom, opop2$pid)]
  opop2$mfmmmffm <- opop2$mfmmmff[match(opop2$mom, opop2$pid)]
  opop2$mfmmfmmm <- opop2$mfmmfmm[match(opop2$mom, opop2$pid)]
  opop2$mfmmfmfm <- opop2$mfmmfmf[match(opop2$mom, opop2$pid)]
  opop2$mfmmffmm <- opop2$mfmmffm[match(opop2$mom, opop2$pid)]
  opop2$mfmmfffm <- opop2$mfmmfff[match(opop2$mom, opop2$pid)]
  opop2$mfmfmmmm <- opop2$mfmfmmm[match(opop2$mom, opop2$pid)]
  opop2$mfmfmmfm <- opop2$mfmfmmf[match(opop2$mom, opop2$pid)]
  opop2$mfmfmfmm <- opop2$mfmfmfm[match(opop2$mom, opop2$pid)]
  opop2$mfmfmffm <- opop2$mfmfmff[match(opop2$mom, opop2$pid)]
  opop2$mfmffmmm <- opop2$mfmffmm[match(opop2$mom, opop2$pid)]
  opop2$mfmffmfm <- opop2$mfmffmf[match(opop2$mom, opop2$pid)]
  opop2$mfmfffmm <- opop2$mfmfffm[match(opop2$mom, opop2$pid)]
  opop2$mfmffffm <- opop2$mfmffff[match(opop2$mom, opop2$pid)]
  opop2$mffmmmmm <- opop2$mffmmmm[match(opop2$mom, opop2$pid)]
  opop2$mffmmmfm <- opop2$mffmmmf[match(opop2$mom, opop2$pid)]
  opop2$mffmmfmm <- opop2$mffmmfm[match(opop2$mom, opop2$pid)]
  opop2$mffmmffm <- opop2$mffmmff[match(opop2$mom, opop2$pid)]
  opop2$mffmfmmm <- opop2$mffmfmm[match(opop2$mom, opop2$pid)]
  opop2$mffmfmfm <- opop2$mffmfmf[match(opop2$mom, opop2$pid)]
  opop2$mffmffmm <- opop2$mffmffm[match(opop2$mom, opop2$pid)]
  opop2$mffmfffm <- opop2$mffmfff[match(opop2$mom, opop2$pid)]
  opop2$mfffmmmm <- opop2$mfffmmm[match(opop2$mom, opop2$pid)]
  opop2$mfffmmfm <- opop2$mfffmmf[match(opop2$mom, opop2$pid)]
  opop2$mfffmfmm <- opop2$mfffmfm[match(opop2$mom, opop2$pid)]
  opop2$mfffmffm <- opop2$mfffmff[match(opop2$mom, opop2$pid)]
  opop2$mffffmmm <- opop2$mffffmm[match(opop2$mom, opop2$pid)]
  opop2$mffffmfm <- opop2$mffffmf[match(opop2$mom, opop2$pid)]
  opop2$mfffffmm <- opop2$mfffffm[match(opop2$mom, opop2$pid)]
  opop2$mffffffm <- opop2$mffffff[match(opop2$mom, opop2$pid)]
  opop2$mmmmmmmf <- opop2$mmmmmmm[match(opop2$pop, opop2$pid)]
  opop2$mmmmmmff <- opop2$mmmmmmf[match(opop2$pop, opop2$pid)]
  opop2$mmmmmfmf <- opop2$mmmmmfm[match(opop2$pop, opop2$pid)]
  opop2$mmmmmfff <- opop2$mmmmmff[match(opop2$pop, opop2$pid)]
  opop2$mmmmfmmf <- opop2$mmmmfmm[match(opop2$pop, opop2$pid)]
  opop2$mmmmfmff <- opop2$mmmmfmf[match(opop2$pop, opop2$pid)]
  opop2$mmmmffmf <- opop2$mmmmffm[match(opop2$pop, opop2$pid)]
  opop2$mmmmffff <- opop2$mmmmfff[match(opop2$pop, opop2$pid)]
  opop2$mmmfmmmf <- opop2$mmmfmmm[match(opop2$pop, opop2$pid)]
  opop2$mmmfmmff <- opop2$mmmfmmf[match(opop2$pop, opop2$pid)]
  opop2$mmmfmfmf <- opop2$mmmfmfm[match(opop2$pop, opop2$pid)]
  opop2$mmmfmfff <- opop2$mmmfmff[match(opop2$pop, opop2$pid)]
  opop2$mmmffmmf <- opop2$mmmffmm[match(opop2$pop, opop2$pid)]
  opop2$mmmffmff <- opop2$mmmffmf[match(opop2$pop, opop2$pid)]
  opop2$mmmfffmf <- opop2$mmmfffm[match(opop2$pop, opop2$pid)]
  opop2$mmmfffff <- opop2$mmmffff[match(opop2$pop, opop2$pid)]
  opop2$mmfmmmmf <- opop2$mmfmmmm[match(opop2$pop, opop2$pid)]
  opop2$mmfmmmff <- opop2$mmfmmmf[match(opop2$pop, opop2$pid)]
  opop2$mmfmmfmf <- opop2$mmfmmfm[match(opop2$pop, opop2$pid)]
  opop2$mmfmmfff <- opop2$mmfmmff[match(opop2$pop, opop2$pid)]
  opop2$mmfmfmmf <- opop2$mmfmfmm[match(opop2$pop, opop2$pid)]
  opop2$mmfmfmff <- opop2$mmfmfmf[match(opop2$pop, opop2$pid)]
  opop2$mmfmffmf <- opop2$mmfmffm[match(opop2$pop, opop2$pid)]
  opop2$mmfmffff <- opop2$mmfmfff[match(opop2$pop, opop2$pid)]
  opop2$mmffmmmf <- opop2$mmffmmm[match(opop2$pop, opop2$pid)]
  opop2$mmffmmff <- opop2$mmffmmf[match(opop2$pop, opop2$pid)]
  opop2$mmffmfmf <- opop2$mmffmfm[match(opop2$pop, opop2$pid)]
  opop2$mmffmfff <- opop2$mmffmff[match(opop2$pop, opop2$pid)]
  opop2$mmfffmmf <- opop2$mmfffmm[match(opop2$pop, opop2$pid)]
  opop2$mmfffmff <- opop2$mmfffmf[match(opop2$pop, opop2$pid)]
  opop2$mmffffmf <- opop2$mmffffm[match(opop2$pop, opop2$pid)]
  opop2$mmffffff <- opop2$mmfffff[match(opop2$pop, opop2$pid)]
  opop2$mfmmmmmf <- opop2$mfmmmmm[match(opop2$pop, opop2$pid)]
  opop2$mfmmmmff <- opop2$mfmmmmf[match(opop2$pop, opop2$pid)]
  opop2$mfmmmfmf <- opop2$mfmmmfm[match(opop2$pop, opop2$pid)]
  opop2$mfmmmfff <- opop2$mfmmmff[match(opop2$pop, opop2$pid)]
  opop2$mfmmfmmf <- opop2$mfmmfmm[match(opop2$pop, opop2$pid)]
  opop2$mfmmfmff <- opop2$mfmmfmf[match(opop2$pop, opop2$pid)]
  opop2$mfmmffmf <- opop2$mfmmffm[match(opop2$pop, opop2$pid)]
  opop2$mfmmffff <- opop2$mfmmfff[match(opop2$pop, opop2$pid)]
  opop2$mfmfmmmf <- opop2$mfmfmmm[match(opop2$pop, opop2$pid)]
  opop2$mfmfmmff <- opop2$mfmfmmf[match(opop2$pop, opop2$pid)]
  opop2$mfmfmfmf <- opop2$mfmfmfm[match(opop2$pop, opop2$pid)]
  opop2$mfmfmfff <- opop2$mfmfmff[match(opop2$pop, opop2$pid)]
  opop2$mfmffmmf <- opop2$mfmffmm[match(opop2$pop, opop2$pid)]
  opop2$mfmffmff <- opop2$mfmffmf[match(opop2$pop, opop2$pid)]
  opop2$mfmfffmf <- opop2$mfmfffm[match(opop2$pop, opop2$pid)]
  opop2$mfmfffff <- opop2$mfmffff[match(opop2$pop, opop2$pid)]
  opop2$mffmmmmf <- opop2$mffmmmm[match(opop2$pop, opop2$pid)]
  opop2$mffmmmff <- opop2$mffmmmf[match(opop2$pop, opop2$pid)]
  opop2$mffmmfmf <- opop2$mffmmfm[match(opop2$pop, opop2$pid)]
  opop2$mffmmfff <- opop2$mffmmff[match(opop2$pop, opop2$pid)]
  opop2$mffmfmmf <- opop2$mffmfmm[match(opop2$pop, opop2$pid)]
  opop2$mffmfmff <- opop2$mffmfmf[match(opop2$pop, opop2$pid)]
  opop2$mffmffmf <- opop2$mffmffm[match(opop2$pop, opop2$pid)]
  opop2$mffmffff <- opop2$mffmfff[match(opop2$pop, opop2$pid)]
  opop2$mfffmmmf <- opop2$mfffmmm[match(opop2$pop, opop2$pid)]
  opop2$mfffmmff <- opop2$mfffmmf[match(opop2$pop, opop2$pid)]
  opop2$mfffmfmf <- opop2$mfffmfm[match(opop2$pop, opop2$pid)]
  opop2$mfffmfff <- opop2$mfffmff[match(opop2$pop, opop2$pid)]
  opop2$mffffmmf <- opop2$mffffmm[match(opop2$pop, opop2$pid)]
  opop2$mffffmff <- opop2$mffffmf[match(opop2$pop, opop2$pid)]
  opop2$mfffffmf <- opop2$mfffffm[match(opop2$pop, opop2$pid)]
  opop2$mfffffff <- opop2$mffffff[match(opop2$pop, opop2$pid)]

  ### 9. Great-great-great-great-great-great-grandfathers
  opop2$fmmmmmmm <- opop2$fmmmmmm[match(opop2$mom, opop2$pid)]
  opop2$fmmmmmfm <- opop2$fmmmmmf[match(opop2$mom, opop2$pid)]
  opop2$fmmmmfmm <- opop2$fmmmmfm[match(opop2$mom, opop2$pid)]
  opop2$fmmmmffm <- opop2$fmmmmff[match(opop2$mom, opop2$pid)]
  opop2$fmmmfmmm <- opop2$fmmmfmm[match(opop2$mom, opop2$pid)]
  opop2$fmmmfmfm <- opop2$fmmmfmf[match(opop2$mom, opop2$pid)]
  opop2$fmmmffmm <- opop2$fmmmffm[match(opop2$mom, opop2$pid)]
  opop2$fmmmfffm <- opop2$fmmmfff[match(opop2$mom, opop2$pid)]
  opop2$fmmfmmmm <- opop2$fmmfmmm[match(opop2$mom, opop2$pid)]
  opop2$fmmfmmfm <- opop2$fmmfmmf[match(opop2$mom, opop2$pid)]
  opop2$fmmfmfmm <- opop2$fmmfmfm[match(opop2$mom, opop2$pid)]
  opop2$fmmfmffm <- opop2$fmmfmff[match(opop2$mom, opop2$pid)]
  opop2$fmmffmmm <- opop2$fmmffmm[match(opop2$mom, opop2$pid)]
  opop2$fmmffmfm <- opop2$fmmffmf[match(opop2$mom, opop2$pid)]
  opop2$fmmfffmm <- opop2$fmmfffm[match(opop2$mom, opop2$pid)]
  opop2$fmmffffm <- opop2$fmmffff[match(opop2$mom, opop2$pid)]
  opop2$fmfmmmmm <- opop2$fmfmmmm[match(opop2$mom, opop2$pid)]
  opop2$fmfmmmfm <- opop2$fmfmmmf[match(opop2$mom, opop2$pid)]
  opop2$fmfmmfmm <- opop2$fmfmmfm[match(opop2$mom, opop2$pid)]
  opop2$fmfmmffm <- opop2$fmfmmff[match(opop2$mom, opop2$pid)]
  opop2$fmfmfmmm <- opop2$fmfmfmm[match(opop2$mom, opop2$pid)]
  opop2$fmfmfmfm <- opop2$fmfmfmf[match(opop2$mom, opop2$pid)]
  opop2$fmfmffmm <- opop2$fmfmffm[match(opop2$mom, opop2$pid)]
  opop2$fmfmfffm <- opop2$fmfmfff[match(opop2$mom, opop2$pid)]
  opop2$fmffmmmm <- opop2$fmffmmm[match(opop2$mom, opop2$pid)]
  opop2$fmffmmfm <- opop2$fmffmmf[match(opop2$mom, opop2$pid)]
  opop2$fmffmfmm <- opop2$fmffmfm[match(opop2$mom, opop2$pid)]
  opop2$fmffmffm <- opop2$fmffmff[match(opop2$mom, opop2$pid)]
  opop2$fmfffmmm <- opop2$fmfffmm[match(opop2$mom, opop2$pid)]
  opop2$fmfffmfm <- opop2$fmfffmf[match(opop2$mom, opop2$pid)]
  opop2$fmffffmm <- opop2$fmffffm[match(opop2$mom, opop2$pid)]
  opop2$fmfffffm <- opop2$fmfffff[match(opop2$mom, opop2$pid)]
  opop2$ffmmmmmm <- opop2$ffmmmmm[match(opop2$mom, opop2$pid)]
  opop2$ffmmmmfm <- opop2$ffmmmmf[match(opop2$mom, opop2$pid)]
  opop2$ffmmmfmm <- opop2$ffmmmfm[match(opop2$mom, opop2$pid)]
  opop2$ffmmmffm <- opop2$ffmmmff[match(opop2$mom, opop2$pid)]
  opop2$ffmmfmmm <- opop2$ffmmfmm[match(opop2$mom, opop2$pid)]
  opop2$ffmmfmfm <- opop2$ffmmfmf[match(opop2$mom, opop2$pid)]
  opop2$ffmmffmm <- opop2$ffmmffm[match(opop2$mom, opop2$pid)]
  opop2$ffmmfffm <- opop2$ffmmfff[match(opop2$mom, opop2$pid)]
  opop2$ffmfmmmm <- opop2$ffmfmmm[match(opop2$mom, opop2$pid)]
  opop2$ffmfmmfm <- opop2$ffmfmmf[match(opop2$mom, opop2$pid)]
  opop2$ffmfmfmm <- opop2$ffmfmfm[match(opop2$mom, opop2$pid)]
  opop2$ffmfmffm <- opop2$ffmfmff[match(opop2$mom, opop2$pid)]
  opop2$ffmffmmm <- opop2$ffmffmm[match(opop2$mom, opop2$pid)]
  opop2$ffmffmfm <- opop2$ffmffmf[match(opop2$mom, opop2$pid)]
  opop2$ffmfffmm <- opop2$ffmfffm[match(opop2$mom, opop2$pid)]
  opop2$ffmffffm <- opop2$ffmffff[match(opop2$mom, opop2$pid)]
  opop2$fffmmmmm <- opop2$fffmmmm[match(opop2$mom, opop2$pid)]
  opop2$fffmmmfm <- opop2$fffmmmf[match(opop2$mom, opop2$pid)]
  opop2$fffmmfmm <- opop2$fffmmfm[match(opop2$mom, opop2$pid)]
  opop2$fffmmffm <- opop2$fffmmff[match(opop2$mom, opop2$pid)]
  opop2$fffmfmmm <- opop2$fffmfmm[match(opop2$mom, opop2$pid)]
  opop2$fffmfmfm <- opop2$fffmfmf[match(opop2$mom, opop2$pid)]
  opop2$fffmffmm <- opop2$fffmffm[match(opop2$mom, opop2$pid)]
  opop2$fffmfffm <- opop2$fffmfff[match(opop2$mom, opop2$pid)]
  opop2$ffffmmmm <- opop2$ffffmmm[match(opop2$mom, opop2$pid)]
  opop2$ffffmmfm <- opop2$ffffmmf[match(opop2$mom, opop2$pid)]
  opop2$ffffmfmm <- opop2$ffffmfm[match(opop2$mom, opop2$pid)]
  opop2$ffffmffm <- opop2$ffffmff[match(opop2$mom, opop2$pid)]
  opop2$fffffmmm <- opop2$fffffmm[match(opop2$mom, opop2$pid)]
  opop2$fffffmfm <- opop2$fffffmf[match(opop2$mom, opop2$pid)]
  opop2$ffffffmm <- opop2$ffffffm[match(opop2$mom, opop2$pid)]
  opop2$fffffffm <- opop2$fffffff[match(opop2$mom, opop2$pid)]
  opop2$fmmmmmmf <- opop2$fmmmmmm[match(opop2$pop, opop2$pid)]
  opop2$fmmmmmff <- opop2$fmmmmmf[match(opop2$pop, opop2$pid)]
  opop2$fmmmmfmf <- opop2$fmmmmfm[match(opop2$pop, opop2$pid)]
  opop2$fmmmmfff <- opop2$fmmmmff[match(opop2$pop, opop2$pid)]
  opop2$fmmmfmmf <- opop2$fmmmfmm[match(opop2$pop, opop2$pid)]
  opop2$fmmmfmff <- opop2$fmmmfmf[match(opop2$pop, opop2$pid)]
  opop2$fmmmffmf <- opop2$fmmmffm[match(opop2$pop, opop2$pid)]
  opop2$fmmmffff <- opop2$fmmmfff[match(opop2$pop, opop2$pid)]
  opop2$fmmfmmmf <- opop2$fmmfmmm[match(opop2$pop, opop2$pid)]
  opop2$fmmfmmff <- opop2$fmmfmmf[match(opop2$pop, opop2$pid)]
  opop2$fmmfmfmf <- opop2$fmmfmfm[match(opop2$pop, opop2$pid)]
  opop2$fmmfmfff <- opop2$fmmfmff[match(opop2$pop, opop2$pid)]
  opop2$fmmffmmf <- opop2$fmmffmm[match(opop2$pop, opop2$pid)]
  opop2$fmmffmff <- opop2$fmmffmf[match(opop2$pop, opop2$pid)]
  opop2$fmmfffmf <- opop2$fmmfffm[match(opop2$pop, opop2$pid)]
  opop2$fmmfffff <- opop2$fmmffff[match(opop2$pop, opop2$pid)]
  opop2$fmfmmmmf <- opop2$fmfmmmm[match(opop2$pop, opop2$pid)]
  opop2$fmfmmmff <- opop2$fmfmmmf[match(opop2$pop, opop2$pid)]
  opop2$fmfmmfmf <- opop2$fmfmmfm[match(opop2$pop, opop2$pid)]
  opop2$fmfmmfff <- opop2$fmfmmff[match(opop2$pop, opop2$pid)]
  opop2$fmfmfmmf <- opop2$fmfmfmm[match(opop2$pop, opop2$pid)]
  opop2$fmfmfmff <- opop2$fmfmfmf[match(opop2$pop, opop2$pid)]
  opop2$fmfmffmf <- opop2$fmfmffm[match(opop2$pop, opop2$pid)]
  opop2$fmfmffff <- opop2$fmfmfff[match(opop2$pop, opop2$pid)]
  opop2$fmffmmmf <- opop2$fmffmmm[match(opop2$pop, opop2$pid)]
  opop2$fmffmmff <- opop2$fmffmmf[match(opop2$pop, opop2$pid)]
  opop2$fmffmfmf <- opop2$fmffmfm[match(opop2$pop, opop2$pid)]
  opop2$fmffmfff <- opop2$fmffmff[match(opop2$pop, opop2$pid)]
  opop2$fmfffmmf <- opop2$fmfffmm[match(opop2$pop, opop2$pid)]
  opop2$fmfffmff <- opop2$fmfffmf[match(opop2$pop, opop2$pid)]
  opop2$fmffffmf <- opop2$fmffffm[match(opop2$pop, opop2$pid)]
  opop2$fmffffff <- opop2$fmfffff[match(opop2$pop, opop2$pid)]
  opop2$ffmmmmmf <- opop2$ffmmmmm[match(opop2$pop, opop2$pid)]
  opop2$ffmmmmff <- opop2$ffmmmmf[match(opop2$pop, opop2$pid)]
  opop2$ffmmmfmf <- opop2$ffmmmfm[match(opop2$pop, opop2$pid)]
  opop2$ffmmmfff <- opop2$ffmmmff[match(opop2$pop, opop2$pid)]
  opop2$ffmmfmmf <- opop2$ffmmfmm[match(opop2$pop, opop2$pid)]
  opop2$ffmmfmff <- opop2$ffmmfmf[match(opop2$pop, opop2$pid)]
  opop2$ffmmffmf <- opop2$ffmmffm[match(opop2$pop, opop2$pid)]
  opop2$ffmmffff <- opop2$ffmmfff[match(opop2$pop, opop2$pid)]
  opop2$ffmfmmmf <- opop2$ffmfmmm[match(opop2$pop, opop2$pid)]
  opop2$ffmfmmff <- opop2$ffmfmmf[match(opop2$pop, opop2$pid)]
  opop2$ffmfmfmf <- opop2$ffmfmfm[match(opop2$pop, opop2$pid)]
  opop2$ffmfmfff <- opop2$ffmfmff[match(opop2$pop, opop2$pid)]
  opop2$ffmffmmf <- opop2$ffmffmm[match(opop2$pop, opop2$pid)]
  opop2$ffmffmff <- opop2$ffmffmf[match(opop2$pop, opop2$pid)]
  opop2$ffmfffmf <- opop2$ffmfffm[match(opop2$pop, opop2$pid)]
  opop2$ffmfffff <- opop2$ffmffff[match(opop2$pop, opop2$pid)]
  opop2$fffmmmmf <- opop2$fffmmmm[match(opop2$pop, opop2$pid)]
  opop2$fffmmmff <- opop2$fffmmmf[match(opop2$pop, opop2$pid)]
  opop2$fffmmfmf <- opop2$fffmmfm[match(opop2$pop, opop2$pid)]
  opop2$fffmmfff <- opop2$fffmmff[match(opop2$pop, opop2$pid)]
  opop2$fffmfmmf <- opop2$fffmfmm[match(opop2$pop, opop2$pid)]
  opop2$fffmfmff <- opop2$fffmfmf[match(opop2$pop, opop2$pid)]
  opop2$fffmffmf <- opop2$fffmffm[match(opop2$pop, opop2$pid)]
  opop2$fffmffff <- opop2$fffmfff[match(opop2$pop, opop2$pid)]
  opop2$ffffmmmf <- opop2$ffffmmm[match(opop2$pop, opop2$pid)]
  opop2$ffffmmff <- opop2$ffffmmf[match(opop2$pop, opop2$pid)]
  opop2$ffffmfmf <- opop2$ffffmfm[match(opop2$pop, opop2$pid)]
  opop2$ffffmfff <- opop2$ffffmff[match(opop2$pop, opop2$pid)]
  opop2$fffffmmf <- opop2$fffffmm[match(opop2$pop, opop2$pid)]
  opop2$fffffmff <- opop2$fffffmf[match(opop2$pop, opop2$pid)]
  opop2$ffffffmf <- opop2$ffffffm[match(opop2$pop, opop2$pid)]
  opop2$ffffffff <- opop2$fffffff[match(opop2$pop, opop2$pid)]
  
  opop2 <- opop2 %>% 
    filter(pid %in% egos) %>% 
    rename(ego = pid) %>% 
    mutate(ego_id = ego) %>% 
    pivot_longer(-ego_id, names_to = "kin_type", values_to = "pid") %>% 
    filter(!is.na(pid)) %>% 
    left_join(select(opop, c(pid, fem, dob, dod, mom, lborn, marid, mstat)), by = "pid") 
  
  return(opop2)
}


## Re-code the kin_type for the ancestors

recode_ancestors <- function(opop_ancestors) {
    opop_ancestors <- opop_ancestors %>% 
    mutate(kin_type = case_when(kin_type %in% c("ego") ~ "ego", 
                                kin_type %in% c("mom", "pop") ~ "parents", 
                                kin_type %in% c("mm", "mf", "fm","ff") ~ "gparents", 
                                kin_type %in% c("mmm", "mfm", "mmf", "mff", 
                                                "fmm", "ffm", "fmf", "fff") ~ "ggparents", 
                                kin_type %in% c("mmmm", "mmfm", "mfmm", "mffm", "mmmf", "mmff", "mfmf", "mfff",
                                                "fmmm", "fmfm", "ffmm", "fffm", "fmmf", "fmff", "ffmf", "ffff") ~ "gggparents", 
                                kin_type %in% c("mmmmm", "mmmfm", "mmfmm", "mmffm", "mfmmm", "mfmfm", "mffmm","mfffm", 
                                                "mmmmf", "mmmff", "mmfmf", "mmfff", "mfmmf", "mfmff", "mffmf", "mffff",
                                                "fmmmm", "fmmfm", "fmfmm", "fmffm", "ffmmm", "ffmfm", "fffmm", "ffffm", 
                                                "fmmmf", "fmmff", "fmfmf", "fmfff", "ffmmf", "ffmff", "fffmf", "fffff") ~ "ggggparents", 
                                kin_type %in% c("mmmmmm", "mmmmfm", "mmmfmm", "mmmffm", "mmfmmm", "mmfmfm", "mmffmm", "mmfffm",
                                                "mfmmmm", "mfmmfm", "mfmfmm", "mfmffm", "mffmmm", "mffmfm", "mfffmm", "mffffm",
                                                "mmmmmf", "mmmmff", "mmmfmf", "mmmfff", "mmfmmf", "mmfmff", "mmffmf", "mmffff",
                                                "mfmmmf", "mfmmff", "mfmfmf", "mfmfff", "mffmmf", "mffmff", "mfffmf", "mfffff",
                                                "fmmmmm", "fmmmfm", "fmmfmm", "fmmffm", "fmfmmm", "fmfmfm", "fmffmm", "fmfffm",
                                                "ffmmmm", "ffmmfm", "ffmfmm", "ffmffm", "fffmmm", "fffmfm", "ffffmm", "fffffm",
                                                "fmmmmf", "fmmmff", "fmmfmf", "fmmfff", "fmfmmf", "fmfmff", "fmffmf", "fmffff",
                                                "ffmmmf", "ffmmff", "ffmfmf", "ffmfff", "fffmmf", "fffmff", "ffffmf", "ffffff") ~ "gggggparents", 
                                kin_type %in% c("mmmmmmm", "mmmmmfm", "mmmmfmm", "mmmmffm", "mmmfmmm", "mmmfmfm", "mmmffmm", "mmmfffm",
                                                "mmfmmmm", "mmfmmfm", "mmfmfmm", "mmfmffm", "mmffmmm", "mmffmfm", "mmfffmm", "mmffffm",
                                                "mfmmmmm", "mfmmmfm", "mfmmfmm", "mfmmffm", "mfmfmmm", "mfmfmfm", "mfmffmm", "mfmfffm",
                                                "mffmmmm", "mffmmfm", "mffmfmm", "mffmffm", "mfffmmm", "mfffmfm", "mffffmm", "mfffffm",
                                                "mmmmmmf", "mmmmmff", "mmmmfmf", "mmmmfff", "mmmfmmf", "mmmfmff", "mmmffmf", "mmmffff",
                                                "mmfmmmf", "mmfmmff", "mmfmfmf", "mmfmfff", "mmffmmf", "mmffmff", "mmfffmf", "mmfffff",
                                                "mfmmmmf", "mfmmmff", "mfmmfmf", "mfmmfff", "mfmfmmf", "mfmfmff", "mfmffmf", "mfmffff",
                                                "mffmmmf", "mffmmff", "mffmfmf", "mffmfff", "mfffmmf", "mfffmff", "mffffmf", "mffffff",
                                                "fmmmmmm", "fmmmmfm", "fmmmfmm", "fmmmffm", "fmmfmmm", "fmmfmfm", "fmmffmm", "fmmfffm",
                                                "fmfmmmm", "fmfmmfm", "fmfmfmm", "fmfmffm", "fmffmmm", "fmffmfm", "fmfffmm", "fmffffm",
                                                "ffmmmmm", "ffmmmfm", "ffmmfmm", "ffmmffm", "ffmfmmm", "ffmfmfm", "ffmffmm", "ffmfffm",
                                                "fffmmmm", "fffmmfm", "fffmfmm", "fffmffm", "ffffmmm", "ffffmfm", "fffffmm", "ffffffm",
                                                "fmmmmmf", "fmmmmff", "fmmmfmf", "fmmmfff", "fmmfmmf", "fmmfmff", "fmmffmf", "fmmffff",
                                                "fmfmmmf", "fmfmmff", "fmfmfmf", "fmfmfff", "fmffmmf", "fmffmff", "fmfffmf", "fmfffff",
                                                "ffmmmmf", "ffmmmff", "ffmmfmf", "ffmmfff", "ffmfmmf", "ffmfmff", "ffmffmf", "ffmffff",
                                                "fffmmmf", "fffmmff", "fffmfmf", "fffmfff", "ffffmmf", "ffffmff", "fffffmf", "fffffff") ~ "ggggggparents", 
                                kin_type %in% c("mmmmmmmm", "mmmmmmfm", "mmmmmfmm", "mmmmmffm", "mmmmfmmm", "mmmmfmfm", "mmmmffmm", "mmmmfffm", 
                                                "mmmfmmmm", "mmmfmmfm", "mmmfmfmm", "mmmfmffm", "mmmffmmm", "mmmffmfm", "mmmfffmm", "mmmffffm", 
                                                "mmfmmmmm", "mmfmmmfm", "mmfmmfmm", "mmfmmffm", "mmfmfmmm", "mmfmfmfm", "mmfmffmm", "mmfmfffm",
                                                "mmffmmmm", "mmffmmfm", "mmffmfmm", "mmffmffm", "mmfffmmm", "mmfffmfm", "mmffffmm", "mmfffffm", 
                                                "mfmmmmmm", "mfmmmmfm", "mfmmmfmm", "mfmmmffm", "mfmmfmmm", "mfmmfmfm", "mfmmffmm", "mfmmfffm", 
                                                "mfmfmmmm", "mfmfmmfm", "mfmfmfmm", "mfmfmffm", "mfmffmmm", "mfmffmfm", "mfmfffmm", "mfmffffm",
                                                "mffmmmmm", "mffmmmfm", "mffmmfmm", "mffmmffm", "mffmfmmm", "mffmfmfm", "mffmffmm", "mffmfffm", 
                                                "mfffmmmm", "mfffmmfm", "mfffmfmm", "mfffmffm", "mffffmmm", "mffffmfm", "mfffffmm", "mffffffm", 
                                                "mmmmmmmf", "mmmmmmff", "mmmmmfmf", "mmmmmfff", "mmmmfmmf", "mmmmfmff", "mmmmffmf", "mmmmffff",
                                                "mmmfmmmf", "mmmfmmff", "mmmfmfmf", "mmmfmfff", "mmmffmmf", "mmmffmff", "mmmfffmf", "mmmfffff", 
                                                "mmfmmmmf", "mmfmmmff", "mmfmmfmf", "mmfmmfff", "mmfmfmmf", "mmfmfmff", "mmfmffmf", "mmfmffff", 
                                                "mmffmmmf", "mmffmmff", "mmffmfmf", "mmffmfff", "mmfffmmf", "mmfffmff", "mmffffmf", "mmffffff",
                                                "mfmmmmmf", "mfmmmmff", "mfmmmfmf", "mfmmmfff", "mfmmfmmf", "mfmmfmff", "mfmmffmf", "mfmmffff", 
                                                "mfmfmmmf", "mfmfmmff", "mfmfmfmf", "mfmfmfff", "mfmffmmf", "mfmffmff", "mfmfffmf", "mfmfffff", 
                                                "mffmmmmf", "mffmmmff", "mffmmfmf", "mffmmfff", "mffmfmmf", "mffmfmff", "mffmffmf", "mffmffff",
                                                "mfffmmmf", "mfffmmff", "mfffmfmf", "mfffmfff", "mffffmmf", "mffffmff", "mfffffmf", "mfffffff", 
                                                "fmmmmmmm", "fmmmmmfm", "fmmmmfmm", "fmmmmffm", "fmmmfmmm", "fmmmfmfm", "fmmmffmm", "fmmmfffm", 
                                                "fmmfmmmm", "fmmfmmfm", "fmmfmfmm", "fmmfmffm", "fmmffmmm", "fmmffmfm", "fmmfffmm", "fmmffffm",
                                                "fmfmmmmm", "fmfmmmfm", "fmfmmfmm", "fmfmmffm", "fmfmfmmm", "fmfmfmfm", "fmfmffmm", "fmfmfffm", 
                                                "fmffmmmm", "fmffmmfm", "fmffmfmm", "fmffmffm", "fmfffmmm", "fmfffmfm", "fmffffmm", "fmfffffm", 
                                                "ffmmmmmm", "ffmmmmfm", "ffmmmfmm", "ffmmmffm", "ffmmfmmm", "ffmmfmfm", "ffmmffmm", "ffmmfffm",
                                                "ffmfmmmm", "ffmfmmfm", "ffmfmfmm", "ffmfmffm", "ffmffmmm", "ffmffmfm", "ffmfffmm", "ffmffffm", 
                                                "fffmmmmm", "fffmmmfm", "fffmmfmm", "fffmmffm", "fffmfmmm", "fffmfmfm", "fffmffmm", "fffmfffm", 
                                                "ffffmmmm", "ffffmmfm", "ffffmfmm", "ffffmffm", "fffffmmm", "fffffmfm", "ffffffmm", "fffffffm",
                                                "fmmmmmmf", "fmmmmmff", "fmmmmfmf", "fmmmmfff", "fmmmfmmf", "fmmmfmff", "fmmmffmf", "fmmmffff", 
                                                "fmmfmmmf", "fmmfmmff", "fmmfmfmf", "fmmfmfff", "fmmffmmf", "fmmffmff", "fmmfffmf", "fmmfffff", 
                                                "fmfmmmmf", "fmfmmmff", "fmfmmfmf", "fmfmmfff", "fmfmfmmf", "fmfmfmff", "fmfmffmf", "fmfmffff", 
                                                "fmffmmmf", "fmffmmff", "fmffmfmf", "fmffmfff", "fmfffmmf", "fmfffmff", "fmffffmf", "fmffffff", 
                                                "ffmmmmmf", "ffmmmmff", "ffmmmfmf", "ffmmmfff", "ffmmfmmf", "ffmmfmff", "ffmmffmf", "ffmmffff", 
                                                "ffmfmmmf", "ffmfmmff", "ffmfmfmf", "ffmfmfff", "ffmffmmf", "ffmffmff", "ffmfffmf", "ffmfffff", 
                                                "fffmmmmf", "fffmmmff", "fffmmfmf", "fffmmfff", "fffmfmmf", "fffmfmff", "fffmffmf", "fffmffff", 
                                                "ffffmmmf", "ffffmmff", "ffffmfmf", "ffffmfff", "fffffmmf", "fffffmff", "ffffffmf", "ffffffff") ~ "gggggggparents", 
                                TRUE ~ NA_character_))
  return(opop_ancestors)
}
