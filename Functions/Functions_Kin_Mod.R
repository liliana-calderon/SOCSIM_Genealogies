#----------------------------------------------------------------------------------------------------
# Function to retrieve members of a kin network of a given ego from SOCSIM microsimulation outputs
# considering direct ancestors up to 8th degree of consanguinity and their offspring

# To run the functions, .opop and omar files must be set in the GlobalEnv

# This is a modified version of the retrieve_kin function created by Mallika Snyder, included in the rsocsim package. 
# Here we add the direct ancestors from 4th to 8th degree of consanguinity and their offspring
# The name of the vectors for grandparents and great-grand parents is also different, 

# Created on 23-09-2022
# Last modified on 13-08-2024
#----------------------------------------------------------------------------------------------------

retrieve_kin_mod <- function (opop = opop, omar = omar, pid, KidsOf = KidsOf) {
  
  "%ni%" <- Negate("%in%")
  zna <- function(x) { return(ifelse(x == 0, NA, x)) }

  opop$spouse <- ifelse(opop$fem, (omar[zna(opop$marid), "hpid"]), 
                        (omar[zna(opop$marid), "wpid"]))
  
  # 3. Grandparents
  opop$mm <- opop$mom[match(opop$mom, opop$pid)]
  opop$mf <- opop$mom[match(opop$pop, opop$pid)]
  opop$fm <- opop$pop[match(opop$mom, opop$pid)]  
  opop$ff <- opop$pop[match(opop$pop, opop$pid)]
  
  # 4. Great-grandparents
  opop$mmm <- opop$mm[match(opop$mom, opop$pid)]
  opop$mfm <- opop$mf[match(opop$mom, opop$pid)]
  opop$mmf <- opop$mm[match(opop$pop, opop$pid)]
  opop$mff <- opop$mf[match(opop$pop, opop$pid)]
  opop$fmm <- opop$fm[match(opop$mom, opop$pid)]
  opop$ffm <- opop$ff[match(opop$mom, opop$pid)] 
  opop$fmf <- opop$fm[match(opop$pop, opop$pid)]
  opop$fff <- opop$ff[match(opop$pop, opop$pid)]
  
  # 5. Great-great-grandparents
  opop$mmmm <- opop$mmm[match(opop$mom, opop$pid)]
  opop$mmfm <- opop$mmf[match(opop$mom, opop$pid)]
  opop$mfmm <- opop$mfm[match(opop$mom, opop$pid)]
  opop$mffm <- opop$mff[match(opop$mom, opop$pid)]
  opop$mmmf <- opop$mmm[match(opop$pop, opop$pid)]
  opop$mmff <- opop$mmf[match(opop$pop, opop$pid)]
  opop$mfmf <- opop$mfm[match(opop$pop, opop$pid)]
  opop$mfff <- opop$mff[match(opop$pop, opop$pid)]
  opop$fmmm <- opop$fmm[match(opop$mom, opop$pid)]
  opop$fmfm <- opop$fmf[match(opop$mom, opop$pid)]
  opop$ffmm <- opop$ffm[match(opop$mom, opop$pid)]
  opop$fffm <- opop$fff[match(opop$mom, opop$pid)]
  opop$fmmf <- opop$fmm[match(opop$pop, opop$pid)]
  opop$fmff <- opop$fmf[match(opop$pop, opop$pid)]
  opop$ffmf <- opop$ffm[match(opop$pop, opop$pid)]
  opop$ffff <- opop$fff[match(opop$pop, opop$pid)]
  
  # 6. Great-great-great-grandparents
  opop$mmmmm <- opop$mmmm[match(opop$mom, opop$pid)]
  opop$mmmfm <- opop$mmmf[match(opop$mom, opop$pid)]
  opop$mmfmm <- opop$mmfm[match(opop$mom, opop$pid)]
  opop$mmffm <- opop$mmff[match(opop$mom, opop$pid)]
  opop$mfmmm <- opop$mfmm[match(opop$mom, opop$pid)]
  opop$mfmfm <- opop$mfmf[match(opop$mom, opop$pid)]
  opop$mffmm <- opop$mffm[match(opop$mom, opop$pid)]
  opop$mfffm <- opop$mfff[match(opop$mom, opop$pid)]
  opop$mmmmf <- opop$mmmm[match(opop$pop, opop$pid)]
  opop$mmmff <- opop$mmmf[match(opop$pop, opop$pid)]
  opop$mmfmf <- opop$mmfm[match(opop$pop, opop$pid)]
  opop$mmfff <- opop$mmff[match(opop$pop, opop$pid)]
  opop$mfmmf <- opop$mfmm[match(opop$pop, opop$pid)]
  opop$mfmff <- opop$mfmf[match(opop$pop, opop$pid)]
  opop$mffmf <- opop$mffm[match(opop$pop, opop$pid)]
  opop$mffff <- opop$mfff[match(opop$pop, opop$pid)]
  opop$fmmmm <- opop$fmmm[match(opop$mom, opop$pid)]
  opop$fmmfm <- opop$fmmf[match(opop$mom, opop$pid)]
  opop$fmfmm <- opop$fmfm[match(opop$mom, opop$pid)]
  opop$fmffm <- opop$fmff[match(opop$mom, opop$pid)]
  opop$ffmmm <- opop$ffmm[match(opop$mom, opop$pid)]
  opop$ffmfm <- opop$ffmf[match(opop$mom, opop$pid)]
  opop$fffmm <- opop$fffm[match(opop$mom, opop$pid)]
  opop$ffffm <- opop$ffff[match(opop$mom, opop$pid)]
  opop$fmmmf <- opop$fmmm[match(opop$pop, opop$pid)]
  opop$fmmff <- opop$fmmf[match(opop$pop, opop$pid)]
  opop$fmfmf <- opop$fmfm[match(opop$pop, opop$pid)]
  opop$fmfff <- opop$fmff[match(opop$pop, opop$pid)]
  opop$ffmmf <- opop$ffmm[match(opop$pop, opop$pid)]
  opop$ffmff <- opop$ffmf[match(opop$pop, opop$pid)]
  opop$fffmf <- opop$fffm[match(opop$pop, opop$pid)]
  opop$fffff <- opop$ffff[match(opop$pop, opop$pid)]
  
  # 7. Great-great-great-great-grandparents
  opop$mmmmmm <- opop$mmmmm[match(opop$mom, opop$pid)]
  opop$mmmmfm <- opop$mmmmf[match(opop$mom, opop$pid)]
  opop$mmmfmm <- opop$mmmfm[match(opop$mom, opop$pid)]
  opop$mmmffm <- opop$mmmff[match(opop$mom, opop$pid)]
  opop$mmfmmm <- opop$mmfmm[match(opop$mom, opop$pid)]
  opop$mmfmfm <- opop$mmfmf[match(opop$mom, opop$pid)]
  opop$mmffmm <- opop$mmffm[match(opop$mom, opop$pid)]
  opop$mmfffm <- opop$mmfff[match(opop$mom, opop$pid)]
  opop$mfmmmm <- opop$mfmmm[match(opop$mom, opop$pid)]
  opop$mfmmfm <- opop$mfmmf[match(opop$mom, opop$pid)]
  opop$mfmfmm <- opop$mfmfm[match(opop$mom, opop$pid)]
  opop$mfmffm <- opop$mfmff[match(opop$mom, opop$pid)]
  opop$mffmmm <- opop$mffmm[match(opop$mom, opop$pid)]
  opop$mffmfm <- opop$mffmf[match(opop$mom, opop$pid)]
  opop$mfffmm <- opop$mfffm[match(opop$mom, opop$pid)]
  opop$mffffm <- opop$mffff[match(opop$mom, opop$pid)]
  opop$mmmmmf <- opop$mmmmm[match(opop$pop, opop$pid)]
  opop$mmmmff <- opop$mmmmf[match(opop$pop, opop$pid)]
  opop$mmmfmf <- opop$mmmfm[match(opop$pop, opop$pid)]
  opop$mmmfff <- opop$mmmff[match(opop$pop, opop$pid)]
  opop$mmfmmf <- opop$mmfmm[match(opop$pop, opop$pid)]
  opop$mmfmff <- opop$mmfmf[match(opop$pop, opop$pid)]
  opop$mmffmf <- opop$mmffm[match(opop$pop, opop$pid)]
  opop$mmffff <- opop$mmfff[match(opop$pop, opop$pid)]
  opop$mfmmmf <- opop$mfmmm[match(opop$pop, opop$pid)]
  opop$mfmmff <- opop$mfmmf[match(opop$pop, opop$pid)]
  opop$mfmfmf <- opop$mfmfm[match(opop$pop, opop$pid)]
  opop$mfmfff <- opop$mfmff[match(opop$pop, opop$pid)]
  opop$mffmmf <- opop$mffmm[match(opop$pop, opop$pid)]
  opop$mffmff <- opop$mffmf[match(opop$pop, opop$pid)]
  opop$mfffmf <- opop$mfffm[match(opop$pop, opop$pid)]
  opop$mfffff <- opop$mffff[match(opop$pop, opop$pid)]
  opop$fmmmmm <- opop$fmmmm[match(opop$mom, opop$pid)]
  opop$fmmmfm <- opop$fmmmf[match(opop$mom, opop$pid)]
  opop$fmmfmm <- opop$fmmfm[match(opop$mom, opop$pid)]
  opop$fmmffm <- opop$fmmff[match(opop$mom, opop$pid)]
  opop$fmfmmm <- opop$fmfmm[match(opop$mom, opop$pid)]
  opop$fmfmfm <- opop$fmfmf[match(opop$mom, opop$pid)]
  opop$fmffmm <- opop$fmffm[match(opop$mom, opop$pid)]
  opop$fmfffm <- opop$fmfff[match(opop$mom, opop$pid)]
  opop$ffmmmm <- opop$ffmmm[match(opop$mom, opop$pid)]
  opop$ffmmfm <- opop$ffmmf[match(opop$mom, opop$pid)]
  opop$ffmfmm <- opop$ffmfm[match(opop$mom, opop$pid)]
  opop$ffmffm <- opop$ffmff[match(opop$mom, opop$pid)]
  opop$fffmmm <- opop$fffmm[match(opop$mom, opop$pid)]
  opop$fffmfm <- opop$fffmf[match(opop$mom, opop$pid)]
  opop$ffffmm <- opop$ffffm[match(opop$mom, opop$pid)]
  opop$fffffm <- opop$fffff[match(opop$mom, opop$pid)]
  opop$fmmmmf <- opop$fmmmm[match(opop$pop, opop$pid)]
  opop$fmmmff <- opop$fmmmf[match(opop$pop, opop$pid)]
  opop$fmmfmf <- opop$fmmfm[match(opop$pop, opop$pid)]
  opop$fmmfff <- opop$fmmff[match(opop$pop, opop$pid)]
  opop$fmfmmf <- opop$fmfmm[match(opop$pop, opop$pid)]
  opop$fmfmff <- opop$fmfmf[match(opop$pop, opop$pid)]
  opop$fmffmf <- opop$fmffm[match(opop$pop, opop$pid)]
  opop$fmffff <- opop$fmfff[match(opop$pop, opop$pid)]
  opop$ffmmmf <- opop$ffmmm[match(opop$pop, opop$pid)]
  opop$ffmmff <- opop$ffmmf[match(opop$pop, opop$pid)]
  opop$ffmfmf <- opop$ffmfm[match(opop$pop, opop$pid)]
  opop$ffmfff <- opop$ffmff[match(opop$pop, opop$pid)]
  opop$fffmmf <- opop$fffmm[match(opop$pop, opop$pid)]
  opop$fffmff <- opop$fffmf[match(opop$pop, opop$pid)]
  opop$ffffmf <- opop$ffffm[match(opop$pop, opop$pid)]
  opop$ffffff <- opop$fffff[match(opop$pop, opop$pid)]
  
  # 8. Great-great-great-great-great-grandparents
  opop$mmmmmmm <- opop$mmmmmm[match(opop$mom, opop$pid)]
  opop$mmmmmfm <- opop$mmmmmf[match(opop$mom, opop$pid)]
  opop$mmmmfmm <- opop$mmmmfm[match(opop$mom, opop$pid)]
  opop$mmmmffm <- opop$mmmmff[match(opop$mom, opop$pid)]
  opop$mmmfmmm <- opop$mmmfmm[match(opop$mom, opop$pid)]
  opop$mmmfmfm <- opop$mmmfmf[match(opop$mom, opop$pid)]
  opop$mmmffmm <- opop$mmmffm[match(opop$mom, opop$pid)]
  opop$mmmfffm <- opop$mmmfff[match(opop$mom, opop$pid)]
  opop$mmfmmmm <- opop$mmfmmm[match(opop$mom, opop$pid)]
  opop$mmfmmfm <- opop$mmfmmf[match(opop$mom, opop$pid)]
  opop$mmfmfmm <- opop$mmfmfm[match(opop$mom, opop$pid)]
  opop$mmfmffm <- opop$mmfmff[match(opop$mom, opop$pid)]
  opop$mmffmmm <- opop$mmffmm[match(opop$mom, opop$pid)]
  opop$mmffmfm <- opop$mmffmf[match(opop$mom, opop$pid)]
  opop$mmfffmm <- opop$mmfffm[match(opop$mom, opop$pid)]
  opop$mmffffm <- opop$mmffff[match(opop$mom, opop$pid)]
  opop$mfmmmmm <- opop$mfmmmm[match(opop$mom, opop$pid)]
  opop$mfmmmfm <- opop$mfmmmf[match(opop$mom, opop$pid)]
  opop$mfmmfmm <- opop$mfmmfm[match(opop$mom, opop$pid)]
  opop$mfmmffm <- opop$mfmmff[match(opop$mom, opop$pid)]
  opop$mfmfmmm <- opop$mfmfmm[match(opop$mom, opop$pid)]
  opop$mfmfmfm <- opop$mfmfmf[match(opop$mom, opop$pid)]
  opop$mfmffmm <- opop$mfmffm[match(opop$mom, opop$pid)]
  opop$mfmfffm <- opop$mfmfff[match(opop$mom, opop$pid)]
  opop$mffmmmm <- opop$mffmmm[match(opop$mom, opop$pid)]
  opop$mffmmfm <- opop$mffmmf[match(opop$mom, opop$pid)]
  opop$mffmfmm <- opop$mffmfm[match(opop$mom, opop$pid)]
  opop$mffmffm <- opop$mffmff[match(opop$mom, opop$pid)]
  opop$mfffmmm <- opop$mfffmm[match(opop$mom, opop$pid)]
  opop$mfffmfm <- opop$mfffmf[match(opop$mom, opop$pid)]
  opop$mffffmm <- opop$mffffm[match(opop$mom, opop$pid)]
  opop$mfffffm <- opop$mfffff[match(opop$mom, opop$pid)]
  opop$mmmmmmf <- opop$mmmmmm[match(opop$pop, opop$pid)]
  opop$mmmmmff <- opop$mmmmmf[match(opop$pop, opop$pid)]
  opop$mmmmfmf <- opop$mmmmfm[match(opop$pop, opop$pid)]
  opop$mmmmfff <- opop$mmmmff[match(opop$pop, opop$pid)]
  opop$mmmfmmf <- opop$mmmfmm[match(opop$pop, opop$pid)]
  opop$mmmfmff <- opop$mmmfmf[match(opop$pop, opop$pid)]
  opop$mmmffmf <- opop$mmmffm[match(opop$pop, opop$pid)]
  opop$mmmffff <- opop$mmmfff[match(opop$pop, opop$pid)]
  opop$mmfmmmf <- opop$mmfmmm[match(opop$pop, opop$pid)]
  opop$mmfmmff <- opop$mmfmmf[match(opop$pop, opop$pid)]
  opop$mmfmfmf <- opop$mmfmfm[match(opop$pop, opop$pid)]
  opop$mmfmfff <- opop$mmfmff[match(opop$pop, opop$pid)]
  opop$mmffmmf <- opop$mmffmm[match(opop$pop, opop$pid)]
  opop$mmffmff <- opop$mmffmf[match(opop$pop, opop$pid)]
  opop$mmfffmf <- opop$mmfffm[match(opop$pop, opop$pid)]
  opop$mmfffff <- opop$mmffff[match(opop$pop, opop$pid)]
  opop$mfmmmmf <- opop$mfmmmm[match(opop$pop, opop$pid)]
  opop$mfmmmff <- opop$mfmmmf[match(opop$pop, opop$pid)]
  opop$mfmmfmf <- opop$mfmmfm[match(opop$pop, opop$pid)]
  opop$mfmmfff <- opop$mfmmff[match(opop$pop, opop$pid)]
  opop$mfmfmmf <- opop$mfmfmm[match(opop$pop, opop$pid)]
  opop$mfmfmff <- opop$mfmfmf[match(opop$pop, opop$pid)]
  opop$mfmffmf <- opop$mfmffm[match(opop$pop, opop$pid)]
  opop$mfmffff <- opop$mfmfff[match(opop$pop, opop$pid)]
  opop$mffmmmf <- opop$mffmmm[match(opop$pop, opop$pid)]
  opop$mffmmff <- opop$mffmmf[match(opop$pop, opop$pid)]
  opop$mffmfmf <- opop$mffmfm[match(opop$pop, opop$pid)]
  opop$mffmfff <- opop$mffmff[match(opop$pop, opop$pid)]
  opop$mfffmmf <- opop$mfffmm[match(opop$pop, opop$pid)]
  opop$mfffmff <- opop$mfffmf[match(opop$pop, opop$pid)]
  opop$mffffmf <- opop$mffffm[match(opop$pop, opop$pid)]
  opop$mffffff <- opop$mfffff[match(opop$pop, opop$pid)]
  opop$fmmmmmm <- opop$fmmmmm[match(opop$mom, opop$pid)]
  opop$fmmmmfm <- opop$fmmmmf[match(opop$mom, opop$pid)]
  opop$fmmmfmm <- opop$fmmmfm[match(opop$mom, opop$pid)]
  opop$fmmmffm <- opop$fmmmff[match(opop$mom, opop$pid)]
  opop$fmmfmmm <- opop$fmmfmm[match(opop$mom, opop$pid)]
  opop$fmmfmfm <- opop$fmmfmf[match(opop$mom, opop$pid)]
  opop$fmmffmm <- opop$fmmffm[match(opop$mom, opop$pid)]
  opop$fmmfffm <- opop$fmmfff[match(opop$mom, opop$pid)]
  opop$fmfmmmm <- opop$fmfmmm[match(opop$mom, opop$pid)]
  opop$fmfmmfm <- opop$fmfmmf[match(opop$mom, opop$pid)]
  opop$fmfmfmm <- opop$fmfmfm[match(opop$mom, opop$pid)]
  opop$fmfmffm <- opop$fmfmff[match(opop$mom, opop$pid)]
  opop$fmffmmm <- opop$fmffmm[match(opop$mom, opop$pid)]
  opop$fmffmfm <- opop$fmffmf[match(opop$mom, opop$pid)]
  opop$fmfffmm <- opop$fmfffm[match(opop$mom, opop$pid)]
  opop$fmffffm <- opop$fmffff[match(opop$mom, opop$pid)]
  opop$ffmmmmm <- opop$ffmmmm[match(opop$mom, opop$pid)]
  opop$ffmmmfm <- opop$ffmmmf[match(opop$mom, opop$pid)]
  opop$ffmmfmm <- opop$ffmmfm[match(opop$mom, opop$pid)]
  opop$ffmmffm <- opop$ffmmff[match(opop$mom, opop$pid)]
  opop$ffmfmmm <- opop$ffmfmm[match(opop$mom, opop$pid)]
  opop$ffmfmfm <- opop$ffmfmf[match(opop$mom, opop$pid)]
  opop$ffmffmm <- opop$ffmffm[match(opop$mom, opop$pid)]
  opop$ffmfffm <- opop$ffmfff[match(opop$mom, opop$pid)]
  opop$fffmmmm <- opop$fffmmm[match(opop$mom, opop$pid)]
  opop$fffmmfm <- opop$fffmmf[match(opop$mom, opop$pid)]
  opop$fffmfmm <- opop$fffmfm[match(opop$mom, opop$pid)]
  opop$fffmffm <- opop$fffmff[match(opop$mom, opop$pid)]
  opop$ffffmmm <- opop$ffffmm[match(opop$mom, opop$pid)]
  opop$ffffmfm <- opop$ffffmf[match(opop$mom, opop$pid)]
  opop$fffffmm <- opop$fffffm[match(opop$mom, opop$pid)]
  opop$ffffffm <- opop$ffffff[match(opop$mom, opop$pid)]
  opop$fmmmmmf <- opop$fmmmmm[match(opop$pop, opop$pid)]
  opop$fmmmmff <- opop$fmmmmf[match(opop$pop, opop$pid)]
  opop$fmmmfmf <- opop$fmmmfm[match(opop$pop, opop$pid)]
  opop$fmmmfff <- opop$fmmmff[match(opop$pop, opop$pid)]
  opop$fmmfmmf <- opop$fmmfmm[match(opop$pop, opop$pid)]
  opop$fmmfmff <- opop$fmmfmf[match(opop$pop, opop$pid)]
  opop$fmmffmf <- opop$fmmffm[match(opop$pop, opop$pid)]
  opop$fmmffff <- opop$fmmfff[match(opop$pop, opop$pid)]
  opop$fmfmmmf <- opop$fmfmmm[match(opop$pop, opop$pid)]
  opop$fmfmmff <- opop$fmfmmf[match(opop$pop, opop$pid)]
  opop$fmfmfmf <- opop$fmfmfm[match(opop$pop, opop$pid)]
  opop$fmfmfff <- opop$fmfmff[match(opop$pop, opop$pid)]
  opop$fmffmmf <- opop$fmffmm[match(opop$pop, opop$pid)]
  opop$fmffmff <- opop$fmffmf[match(opop$pop, opop$pid)]
  opop$fmfffmf <- opop$fmfffm[match(opop$pop, opop$pid)]
  opop$fmfffff <- opop$fmffff[match(opop$pop, opop$pid)]
  opop$ffmmmmf <- opop$ffmmmm[match(opop$pop, opop$pid)]
  opop$ffmmmff <- opop$ffmmmf[match(opop$pop, opop$pid)]
  opop$ffmmfmf <- opop$ffmmfm[match(opop$pop, opop$pid)]
  opop$ffmmfff <- opop$ffmmff[match(opop$pop, opop$pid)]
  opop$ffmfmmf <- opop$ffmfmm[match(opop$pop, opop$pid)]
  opop$ffmfmff <- opop$ffmfmf[match(opop$pop, opop$pid)]
  opop$ffmffmf <- opop$ffmffm[match(opop$pop, opop$pid)]
  opop$ffmffff <- opop$ffmfff[match(opop$pop, opop$pid)]
  opop$fffmmmf <- opop$fffmmm[match(opop$pop, opop$pid)]
  opop$fffmmff <- opop$fffmmf[match(opop$pop, opop$pid)]
  opop$fffmfmf <- opop$fffmfm[match(opop$pop, opop$pid)]
  opop$fffmfff <- opop$fffmff[match(opop$pop, opop$pid)]
  opop$ffffmmf <- opop$ffffmm[match(opop$pop, opop$pid)]
  opop$ffffmff <- opop$ffffmf[match(opop$pop, opop$pid)]
  opop$fffffmf <- opop$fffffm[match(opop$pop, opop$pid)]
  opop$fffffff <- opop$ffffff[match(opop$pop, opop$pid)]
  
  # 9. Great-great-great-great-great-great-grandparents
  opop$mmmmmmmm <- opop$mmmmmmm[match(opop$mom, opop$pid)]
  opop$mmmmmmfm <- opop$mmmmmmf[match(opop$mom, opop$pid)]
  opop$mmmmmfmm <- opop$mmmmmfm[match(opop$mom, opop$pid)]
  opop$mmmmmffm <- opop$mmmmmff[match(opop$mom, opop$pid)]
  opop$mmmmfmmm <- opop$mmmmfmm[match(opop$mom, opop$pid)]
  opop$mmmmfmfm <- opop$mmmmfmf[match(opop$mom, opop$pid)]
  opop$mmmmffmm <- opop$mmmmffm[match(opop$mom, opop$pid)]
  opop$mmmmfffm <- opop$mmmmfff[match(opop$mom, opop$pid)]
  opop$mmmfmmmm <- opop$mmmfmmm[match(opop$mom, opop$pid)]
  opop$mmmfmmfm <- opop$mmmfmmf[match(opop$mom, opop$pid)]
  opop$mmmfmfmm <- opop$mmmfmfm[match(opop$mom, opop$pid)]
  opop$mmmfmffm <- opop$mmmfmff[match(opop$mom, opop$pid)]
  opop$mmmffmmm <- opop$mmmffmm[match(opop$mom, opop$pid)]
  opop$mmmffmfm <- opop$mmmffmf[match(opop$mom, opop$pid)]
  opop$mmmfffmm <- opop$mmmfffm[match(opop$mom, opop$pid)]
  opop$mmmffffm <- opop$mmmffff[match(opop$mom, opop$pid)]
  opop$mmfmmmmm <- opop$mmfmmmm[match(opop$mom, opop$pid)]
  opop$mmfmmmfm <- opop$mmfmmmf[match(opop$mom, opop$pid)]
  opop$mmfmmfmm <- opop$mmfmmfm[match(opop$mom, opop$pid)]
  opop$mmfmmffm <- opop$mmfmmff[match(opop$mom, opop$pid)]
  opop$mmfmfmmm <- opop$mmfmfmm[match(opop$mom, opop$pid)]
  opop$mmfmfmfm <- opop$mmfmfmf[match(opop$mom, opop$pid)]
  opop$mmfmffmm <- opop$mmfmffm[match(opop$mom, opop$pid)]
  opop$mmfmfffm <- opop$mmfmfff[match(opop$mom, opop$pid)]
  opop$mmffmmmm <- opop$mmffmmm[match(opop$mom, opop$pid)]
  opop$mmffmmfm <- opop$mmffmmf[match(opop$mom, opop$pid)]
  opop$mmffmfmm <- opop$mmffmfm[match(opop$mom, opop$pid)]
  opop$mmffmffm <- opop$mmffmff[match(opop$mom, opop$pid)]
  opop$mmfffmmm <- opop$mmfffmm[match(opop$mom, opop$pid)]
  opop$mmfffmfm <- opop$mmfffmf[match(opop$mom, opop$pid)]
  opop$mmffffmm <- opop$mmffffm[match(opop$mom, opop$pid)]
  opop$mmfffffm <- opop$mmfffff[match(opop$mom, opop$pid)]
  opop$mfmmmmmm <- opop$mfmmmmm[match(opop$mom, opop$pid)]
  opop$mfmmmmfm <- opop$mfmmmmf[match(opop$mom, opop$pid)]
  opop$mfmmmfmm <- opop$mfmmmfm[match(opop$mom, opop$pid)]
  opop$mfmmmffm <- opop$mfmmmff[match(opop$mom, opop$pid)]
  opop$mfmmfmmm <- opop$mfmmfmm[match(opop$mom, opop$pid)]
  opop$mfmmfmfm <- opop$mfmmfmf[match(opop$mom, opop$pid)]
  opop$mfmmffmm <- opop$mfmmffm[match(opop$mom, opop$pid)]
  opop$mfmmfffm <- opop$mfmmfff[match(opop$mom, opop$pid)]
  opop$mfmfmmmm <- opop$mfmfmmm[match(opop$mom, opop$pid)]
  opop$mfmfmmfm <- opop$mfmfmmf[match(opop$mom, opop$pid)]
  opop$mfmfmfmm <- opop$mfmfmfm[match(opop$mom, opop$pid)]
  opop$mfmfmffm <- opop$mfmfmff[match(opop$mom, opop$pid)]
  opop$mfmffmmm <- opop$mfmffmm[match(opop$mom, opop$pid)]
  opop$mfmffmfm <- opop$mfmffmf[match(opop$mom, opop$pid)]
  opop$mfmfffmm <- opop$mfmfffm[match(opop$mom, opop$pid)]
  opop$mfmffffm <- opop$mfmffff[match(opop$mom, opop$pid)]
  opop$mffmmmmm <- opop$mffmmmm[match(opop$mom, opop$pid)]
  opop$mffmmmfm <- opop$mffmmmf[match(opop$mom, opop$pid)]
  opop$mffmmfmm <- opop$mffmmfm[match(opop$mom, opop$pid)]
  opop$mffmmffm <- opop$mffmmff[match(opop$mom, opop$pid)]
  opop$mffmfmmm <- opop$mffmfmm[match(opop$mom, opop$pid)]
  opop$mffmfmfm <- opop$mffmfmf[match(opop$mom, opop$pid)]
  opop$mffmffmm <- opop$mffmffm[match(opop$mom, opop$pid)]
  opop$mffmfffm <- opop$mffmfff[match(opop$mom, opop$pid)]
  opop$mfffmmmm <- opop$mfffmmm[match(opop$mom, opop$pid)]
  opop$mfffmmfm <- opop$mfffmmf[match(opop$mom, opop$pid)]
  opop$mfffmfmm <- opop$mfffmfm[match(opop$mom, opop$pid)]
  opop$mfffmffm <- opop$mfffmff[match(opop$mom, opop$pid)]
  opop$mffffmmm <- opop$mffffmm[match(opop$mom, opop$pid)]
  opop$mffffmfm <- opop$mffffmf[match(opop$mom, opop$pid)]
  opop$mfffffmm <- opop$mfffffm[match(opop$mom, opop$pid)]
  opop$mffffffm <- opop$mffffff[match(opop$mom, opop$pid)]
  opop$mmmmmmmf <- opop$mmmmmmm[match(opop$pop, opop$pid)]
  opop$mmmmmmff <- opop$mmmmmmf[match(opop$pop, opop$pid)]
  opop$mmmmmfmf <- opop$mmmmmfm[match(opop$pop, opop$pid)]
  opop$mmmmmfff <- opop$mmmmmff[match(opop$pop, opop$pid)]
  opop$mmmmfmmf <- opop$mmmmfmm[match(opop$pop, opop$pid)]
  opop$mmmmfmff <- opop$mmmmfmf[match(opop$pop, opop$pid)]
  opop$mmmmffmf <- opop$mmmmffm[match(opop$pop, opop$pid)]
  opop$mmmmffff <- opop$mmmmfff[match(opop$pop, opop$pid)]
  opop$mmmfmmmf <- opop$mmmfmmm[match(opop$pop, opop$pid)]
  opop$mmmfmmff <- opop$mmmfmmf[match(opop$pop, opop$pid)]
  opop$mmmfmfmf <- opop$mmmfmfm[match(opop$pop, opop$pid)]
  opop$mmmfmfff <- opop$mmmfmff[match(opop$pop, opop$pid)]
  opop$mmmffmmf <- opop$mmmffmm[match(opop$pop, opop$pid)]
  opop$mmmffmff <- opop$mmmffmf[match(opop$pop, opop$pid)]
  opop$mmmfffmf <- opop$mmmfffm[match(opop$pop, opop$pid)]
  opop$mmmfffff <- opop$mmmffff[match(opop$pop, opop$pid)]
  opop$mmfmmmmf <- opop$mmfmmmm[match(opop$pop, opop$pid)]
  opop$mmfmmmff <- opop$mmfmmmf[match(opop$pop, opop$pid)]
  opop$mmfmmfmf <- opop$mmfmmfm[match(opop$pop, opop$pid)]
  opop$mmfmmfff <- opop$mmfmmff[match(opop$pop, opop$pid)]
  opop$mmfmfmmf <- opop$mmfmfmm[match(opop$pop, opop$pid)]
  opop$mmfmfmff <- opop$mmfmfmf[match(opop$pop, opop$pid)]
  opop$mmfmffmf <- opop$mmfmffm[match(opop$pop, opop$pid)]
  opop$mmfmffff <- opop$mmfmfff[match(opop$pop, opop$pid)]
  opop$mmffmmmf <- opop$mmffmmm[match(opop$pop, opop$pid)]
  opop$mmffmmff <- opop$mmffmmf[match(opop$pop, opop$pid)]
  opop$mmffmfmf <- opop$mmffmfm[match(opop$pop, opop$pid)]
  opop$mmffmfff <- opop$mmffmff[match(opop$pop, opop$pid)]
  opop$mmfffmmf <- opop$mmfffmm[match(opop$pop, opop$pid)]
  opop$mmfffmff <- opop$mmfffmf[match(opop$pop, opop$pid)]
  opop$mmffffmf <- opop$mmffffm[match(opop$pop, opop$pid)]
  opop$mmffffff <- opop$mmfffff[match(opop$pop, opop$pid)]
  opop$mfmmmmmf <- opop$mfmmmmm[match(opop$pop, opop$pid)]
  opop$mfmmmmff <- opop$mfmmmmf[match(opop$pop, opop$pid)]
  opop$mfmmmfmf <- opop$mfmmmfm[match(opop$pop, opop$pid)]
  opop$mfmmmfff <- opop$mfmmmff[match(opop$pop, opop$pid)]
  opop$mfmmfmmf <- opop$mfmmfmm[match(opop$pop, opop$pid)]
  opop$mfmmfmff <- opop$mfmmfmf[match(opop$pop, opop$pid)]
  opop$mfmmffmf <- opop$mfmmffm[match(opop$pop, opop$pid)]
  opop$mfmmffff <- opop$mfmmfff[match(opop$pop, opop$pid)]
  opop$mfmfmmmf <- opop$mfmfmmm[match(opop$pop, opop$pid)]
  opop$mfmfmmff <- opop$mfmfmmf[match(opop$pop, opop$pid)]
  opop$mfmfmfmf <- opop$mfmfmfm[match(opop$pop, opop$pid)]
  opop$mfmfmfff <- opop$mfmfmff[match(opop$pop, opop$pid)]
  opop$mfmffmmf <- opop$mfmffmm[match(opop$pop, opop$pid)]
  opop$mfmffmff <- opop$mfmffmf[match(opop$pop, opop$pid)]
  opop$mfmfffmf <- opop$mfmfffm[match(opop$pop, opop$pid)]
  opop$mfmfffff <- opop$mfmffff[match(opop$pop, opop$pid)]
  opop$mffmmmmf <- opop$mffmmmm[match(opop$pop, opop$pid)]
  opop$mffmmmff <- opop$mffmmmf[match(opop$pop, opop$pid)]
  opop$mffmmfmf <- opop$mffmmfm[match(opop$pop, opop$pid)]
  opop$mffmmfff <- opop$mffmmff[match(opop$pop, opop$pid)]
  opop$mffmfmmf <- opop$mffmfmm[match(opop$pop, opop$pid)]
  opop$mffmfmff <- opop$mffmfmf[match(opop$pop, opop$pid)]
  opop$mffmffmf <- opop$mffmffm[match(opop$pop, opop$pid)]
  opop$mffmffff <- opop$mffmfff[match(opop$pop, opop$pid)]
  opop$mfffmmmf <- opop$mfffmmm[match(opop$pop, opop$pid)]
  opop$mfffmmff <- opop$mfffmmf[match(opop$pop, opop$pid)]
  opop$mfffmfmf <- opop$mfffmfm[match(opop$pop, opop$pid)]
  opop$mfffmfff <- opop$mfffmff[match(opop$pop, opop$pid)]
  opop$mffffmmf <- opop$mffffmm[match(opop$pop, opop$pid)]
  opop$mffffmff <- opop$mffffmf[match(opop$pop, opop$pid)]
  opop$mfffffmf <- opop$mfffffm[match(opop$pop, opop$pid)]
  opop$mfffffff <- opop$mffffff[match(opop$pop, opop$pid)]
  opop$fmmmmmmm <- opop$fmmmmmm[match(opop$mom, opop$pid)]
  opop$fmmmmmfm <- opop$fmmmmmf[match(opop$mom, opop$pid)]
  opop$fmmmmfmm <- opop$fmmmmfm[match(opop$mom, opop$pid)]
  opop$fmmmmffm <- opop$fmmmmff[match(opop$mom, opop$pid)]
  opop$fmmmfmmm <- opop$fmmmfmm[match(opop$mom, opop$pid)]
  opop$fmmmfmfm <- opop$fmmmfmf[match(opop$mom, opop$pid)]
  opop$fmmmffmm <- opop$fmmmffm[match(opop$mom, opop$pid)]
  opop$fmmmfffm <- opop$fmmmfff[match(opop$mom, opop$pid)]
  opop$fmmfmmmm <- opop$fmmfmmm[match(opop$mom, opop$pid)]
  opop$fmmfmmfm <- opop$fmmfmmf[match(opop$mom, opop$pid)]
  opop$fmmfmfmm <- opop$fmmfmfm[match(opop$mom, opop$pid)]
  opop$fmmfmffm <- opop$fmmfmff[match(opop$mom, opop$pid)]
  opop$fmmffmmm <- opop$fmmffmm[match(opop$mom, opop$pid)]
  opop$fmmffmfm <- opop$fmmffmf[match(opop$mom, opop$pid)]
  opop$fmmfffmm <- opop$fmmfffm[match(opop$mom, opop$pid)]
  opop$fmmffffm <- opop$fmmffff[match(opop$mom, opop$pid)]
  opop$fmfmmmmm <- opop$fmfmmmm[match(opop$mom, opop$pid)]
  opop$fmfmmmfm <- opop$fmfmmmf[match(opop$mom, opop$pid)]
  opop$fmfmmfmm <- opop$fmfmmfm[match(opop$mom, opop$pid)]
  opop$fmfmmffm <- opop$fmfmmff[match(opop$mom, opop$pid)]
  opop$fmfmfmmm <- opop$fmfmfmm[match(opop$mom, opop$pid)]
  opop$fmfmfmfm <- opop$fmfmfmf[match(opop$mom, opop$pid)]
  opop$fmfmffmm <- opop$fmfmffm[match(opop$mom, opop$pid)]
  opop$fmfmfffm <- opop$fmfmfff[match(opop$mom, opop$pid)]
  opop$fmffmmmm <- opop$fmffmmm[match(opop$mom, opop$pid)]
  opop$fmffmmfm <- opop$fmffmmf[match(opop$mom, opop$pid)]
  opop$fmffmfmm <- opop$fmffmfm[match(opop$mom, opop$pid)]
  opop$fmffmffm <- opop$fmffmff[match(opop$mom, opop$pid)]
  opop$fmfffmmm <- opop$fmfffmm[match(opop$mom, opop$pid)]
  opop$fmfffmfm <- opop$fmfffmf[match(opop$mom, opop$pid)]
  opop$fmffffmm <- opop$fmffffm[match(opop$mom, opop$pid)]
  opop$fmfffffm <- opop$fmfffff[match(opop$mom, opop$pid)]
  opop$ffmmmmmm <- opop$ffmmmmm[match(opop$mom, opop$pid)]
  opop$ffmmmmfm <- opop$ffmmmmf[match(opop$mom, opop$pid)]
  opop$ffmmmfmm <- opop$ffmmmfm[match(opop$mom, opop$pid)]
  opop$ffmmmffm <- opop$ffmmmff[match(opop$mom, opop$pid)]
  opop$ffmmfmmm <- opop$ffmmfmm[match(opop$mom, opop$pid)]
  opop$ffmmfmfm <- opop$ffmmfmf[match(opop$mom, opop$pid)]
  opop$ffmmffmm <- opop$ffmmffm[match(opop$mom, opop$pid)]
  opop$ffmmfffm <- opop$ffmmfff[match(opop$mom, opop$pid)]
  opop$ffmfmmmm <- opop$ffmfmmm[match(opop$mom, opop$pid)]
  opop$ffmfmmfm <- opop$ffmfmmf[match(opop$mom, opop$pid)]
  opop$ffmfmfmm <- opop$ffmfmfm[match(opop$mom, opop$pid)]
  opop$ffmfmffm <- opop$ffmfmff[match(opop$mom, opop$pid)]
  opop$ffmffmmm <- opop$ffmffmm[match(opop$mom, opop$pid)]
  opop$ffmffmfm <- opop$ffmffmf[match(opop$mom, opop$pid)]
  opop$ffmfffmm <- opop$ffmfffm[match(opop$mom, opop$pid)]
  opop$ffmffffm <- opop$ffmffff[match(opop$mom, opop$pid)]
  opop$fffmmmmm <- opop$fffmmmm[match(opop$mom, opop$pid)]
  opop$fffmmmfm <- opop$fffmmmf[match(opop$mom, opop$pid)]
  opop$fffmmfmm <- opop$fffmmfm[match(opop$mom, opop$pid)]
  opop$fffmmffm <- opop$fffmmff[match(opop$mom, opop$pid)]
  opop$fffmfmmm <- opop$fffmfmm[match(opop$mom, opop$pid)]
  opop$fffmfmfm <- opop$fffmfmf[match(opop$mom, opop$pid)]
  opop$fffmffmm <- opop$fffmffm[match(opop$mom, opop$pid)]
  opop$fffmfffm <- opop$fffmfff[match(opop$mom, opop$pid)]
  opop$ffffmmmm <- opop$ffffmmm[match(opop$mom, opop$pid)]
  opop$ffffmmfm <- opop$ffffmmf[match(opop$mom, opop$pid)]
  opop$ffffmfmm <- opop$ffffmfm[match(opop$mom, opop$pid)]
  opop$ffffmffm <- opop$ffffmff[match(opop$mom, opop$pid)]
  opop$fffffmmm <- opop$fffffmm[match(opop$mom, opop$pid)]
  opop$fffffmfm <- opop$fffffmf[match(opop$mom, opop$pid)]
  opop$ffffffmm <- opop$ffffffm[match(opop$mom, opop$pid)]
  opop$fffffffm <- opop$fffffff[match(opop$mom, opop$pid)]
  opop$fmmmmmmf <- opop$fmmmmmm[match(opop$pop, opop$pid)]
  opop$fmmmmmff <- opop$fmmmmmf[match(opop$pop, opop$pid)]
  opop$fmmmmfmf <- opop$fmmmmfm[match(opop$pop, opop$pid)]
  opop$fmmmmfff <- opop$fmmmmff[match(opop$pop, opop$pid)]
  opop$fmmmfmmf <- opop$fmmmfmm[match(opop$pop, opop$pid)]
  opop$fmmmfmff <- opop$fmmmfmf[match(opop$pop, opop$pid)]
  opop$fmmmffmf <- opop$fmmmffm[match(opop$pop, opop$pid)]
  opop$fmmmffff <- opop$fmmmfff[match(opop$pop, opop$pid)]
  opop$fmmfmmmf <- opop$fmmfmmm[match(opop$pop, opop$pid)]
  opop$fmmfmmff <- opop$fmmfmmf[match(opop$pop, opop$pid)]
  opop$fmmfmfmf <- opop$fmmfmfm[match(opop$pop, opop$pid)]
  opop$fmmfmfff <- opop$fmmfmff[match(opop$pop, opop$pid)]
  opop$fmmffmmf <- opop$fmmffmm[match(opop$pop, opop$pid)]
  opop$fmmffmff <- opop$fmmffmf[match(opop$pop, opop$pid)]
  opop$fmmfffmf <- opop$fmmfffm[match(opop$pop, opop$pid)]
  opop$fmmfffff <- opop$fmmffff[match(opop$pop, opop$pid)]
  opop$fmfmmmmf <- opop$fmfmmmm[match(opop$pop, opop$pid)]
  opop$fmfmmmff <- opop$fmfmmmf[match(opop$pop, opop$pid)]
  opop$fmfmmfmf <- opop$fmfmmfm[match(opop$pop, opop$pid)]
  opop$fmfmmfff <- opop$fmfmmff[match(opop$pop, opop$pid)]
  opop$fmfmfmmf <- opop$fmfmfmm[match(opop$pop, opop$pid)]
  opop$fmfmfmff <- opop$fmfmfmf[match(opop$pop, opop$pid)]
  opop$fmfmffmf <- opop$fmfmffm[match(opop$pop, opop$pid)]
  opop$fmfmffff <- opop$fmfmfff[match(opop$pop, opop$pid)]
  opop$fmffmmmf <- opop$fmffmmm[match(opop$pop, opop$pid)]
  opop$fmffmmff <- opop$fmffmmf[match(opop$pop, opop$pid)]
  opop$fmffmfmf <- opop$fmffmfm[match(opop$pop, opop$pid)]
  opop$fmffmfff <- opop$fmffmff[match(opop$pop, opop$pid)]
  opop$fmfffmmf <- opop$fmfffmm[match(opop$pop, opop$pid)]
  opop$fmfffmff <- opop$fmfffmf[match(opop$pop, opop$pid)]
  opop$fmffffmf <- opop$fmffffm[match(opop$pop, opop$pid)]
  opop$fmffffff <- opop$fmfffff[match(opop$pop, opop$pid)]
  opop$ffmmmmmf <- opop$ffmmmmm[match(opop$pop, opop$pid)]
  opop$ffmmmmff <- opop$ffmmmmf[match(opop$pop, opop$pid)]
  opop$ffmmmfmf <- opop$ffmmmfm[match(opop$pop, opop$pid)]
  opop$ffmmmfff <- opop$ffmmmff[match(opop$pop, opop$pid)]
  opop$ffmmfmmf <- opop$ffmmfmm[match(opop$pop, opop$pid)]
  opop$ffmmfmff <- opop$ffmmfmf[match(opop$pop, opop$pid)]
  opop$ffmmffmf <- opop$ffmmffm[match(opop$pop, opop$pid)]
  opop$ffmmffff <- opop$ffmmfff[match(opop$pop, opop$pid)]
  opop$ffmfmmmf <- opop$ffmfmmm[match(opop$pop, opop$pid)]
  opop$ffmfmmff <- opop$ffmfmmf[match(opop$pop, opop$pid)]
  opop$ffmfmfmf <- opop$ffmfmfm[match(opop$pop, opop$pid)]
  opop$ffmfmfff <- opop$ffmfmff[match(opop$pop, opop$pid)]
  opop$ffmffmmf <- opop$ffmffmm[match(opop$pop, opop$pid)]
  opop$ffmffmff <- opop$ffmffmf[match(opop$pop, opop$pid)]
  opop$ffmfffmf <- opop$ffmfffm[match(opop$pop, opop$pid)]
  opop$ffmfffff <- opop$ffmffff[match(opop$pop, opop$pid)]
  opop$fffmmmmf <- opop$fffmmmm[match(opop$pop, opop$pid)]
  opop$fffmmmff <- opop$fffmmmf[match(opop$pop, opop$pid)]
  opop$fffmmfmf <- opop$fffmmfm[match(opop$pop, opop$pid)]
  opop$fffmmfff <- opop$fffmmff[match(opop$pop, opop$pid)]
  opop$fffmfmmf <- opop$fffmfmm[match(opop$pop, opop$pid)]
  opop$fffmfmff <- opop$fffmfmf[match(opop$pop, opop$pid)]
  opop$fffmffmf <- opop$fffmffm[match(opop$pop, opop$pid)]
  opop$fffmffff <- opop$fffmfff[match(opop$pop, opop$pid)]
  opop$ffffmmmf <- opop$ffffmmm[match(opop$pop, opop$pid)]
  opop$ffffmmff <- opop$ffffmmf[match(opop$pop, opop$pid)]
  opop$ffffmfmf <- opop$ffffmfm[match(opop$pop, opop$pid)]
  opop$ffffmfff <- opop$ffffmff[match(opop$pop, opop$pid)]
  opop$fffffmmf <- opop$fffffmm[match(opop$pop, opop$pid)]
  opop$fffffmff <- opop$fffffmf[match(opop$pop, opop$pid)]
  opop$ffffffmf <- opop$ffffffm[match(opop$pop, opop$pid)]
  opop$ffffffff <- opop$fffffff[match(opop$pop, opop$pid)]
  
  # 10. Great-great-great-great-great-great-great-grandparents
  opop$mmmmmmmmm <- opop$mmmmmmmm[match(opop$mom, opop$pid)]
  opop$mmmmmmfmm <- opop$mmmmmmfm[match(opop$mom, opop$pid)]
  opop$mmmmmfmmm <- opop$mmmmmfmm[match(opop$mom, opop$pid)]
  opop$mmmmmffmm <- opop$mmmmmffm[match(opop$mom, opop$pid)]
  opop$mmmmfmmmm <- opop$mmmmfmmm[match(opop$mom, opop$pid)]
  opop$mmmmfmfmm <- opop$mmmmfmfm[match(opop$mom, opop$pid)]
  opop$mmmmffmmm <- opop$mmmmffmm[match(opop$mom, opop$pid)]
  opop$mmmmfffmm <- opop$mmmmfffm[match(opop$mom, opop$pid)]
  opop$mmmfmmmmm <- opop$mmmfmmmm[match(opop$mom, opop$pid)]
  opop$mmmfmmfmm <- opop$mmmfmmfm[match(opop$mom, opop$pid)]
  opop$mmmfmfmmm <- opop$mmmfmfmm[match(opop$mom, opop$pid)]
  opop$mmmfmffmm <- opop$mmmfmffm[match(opop$mom, opop$pid)]
  opop$mmmffmmmm <- opop$mmmffmmm[match(opop$mom, opop$pid)]
  opop$mmmffmfmm <- opop$mmmffmfm[match(opop$mom, opop$pid)]
  opop$mmmfffmmm <- opop$mmmfffmm[match(opop$mom, opop$pid)]
  opop$mmmffffmm <- opop$mmmffffm[match(opop$mom, opop$pid)]
  opop$mmfmmmmmm <- opop$mmfmmmmm[match(opop$mom, opop$pid)]
  opop$mmfmmmfmm <- opop$mmfmmmfm[match(opop$mom, opop$pid)]
  opop$mmfmmfmmm <- opop$mmfmmfmm[match(opop$mom, opop$pid)]
  opop$mmfmmffmm <- opop$mmfmmffm[match(opop$mom, opop$pid)]
  opop$mmfmfmmmm <- opop$mmfmfmmm[match(opop$mom, opop$pid)]
  opop$mmfmfmfmm <- opop$mmfmfmfm[match(opop$mom, opop$pid)]
  opop$mmfmffmmm <- opop$mmfmffmm[match(opop$mom, opop$pid)]
  opop$mmfmfffmm <- opop$mmfmfffm[match(opop$mom, opop$pid)]
  opop$mmffmmmmm <- opop$mmffmmmm[match(opop$mom, opop$pid)]
  opop$mmffmmfmm <- opop$mmffmmfm[match(opop$mom, opop$pid)]
  opop$mmffmfmmm <- opop$mmffmfmm[match(opop$mom, opop$pid)]
  opop$mmffmffmm <- opop$mmffmffm[match(opop$mom, opop$pid)]
  opop$mmfffmmmm <- opop$mmfffmmm[match(opop$mom, opop$pid)]
  opop$mmfffmfmm <- opop$mmfffmfm[match(opop$mom, opop$pid)]
  opop$mmffffmmm <- opop$mmffffmm[match(opop$mom, opop$pid)]
  opop$mmfffffmm <- opop$mmfffffm[match(opop$mom, opop$pid)]
  opop$mfmmmmmmm <- opop$mfmmmmmm[match(opop$mom, opop$pid)]
  opop$mfmmmmfmm <- opop$mfmmmmfm[match(opop$mom, opop$pid)]
  opop$mfmmmfmmm <- opop$mfmmmfmm[match(opop$mom, opop$pid)]
  opop$mfmmmffmm <- opop$mfmmmffm[match(opop$mom, opop$pid)]
  opop$mfmmfmmmm <- opop$mfmmfmmm[match(opop$mom, opop$pid)]
  opop$mfmmfmfmm <- opop$mfmmfmfm[match(opop$mom, opop$pid)]
  opop$mfmmffmmm <- opop$mfmmffmm[match(opop$mom, opop$pid)]
  opop$mfmmfffmm <- opop$mfmmfffm[match(opop$mom, opop$pid)]
  opop$mfmfmmmmm <- opop$mfmfmmmm[match(opop$mom, opop$pid)]
  opop$mfmfmmfmm <- opop$mfmfmmfm[match(opop$mom, opop$pid)]
  opop$mfmfmfmmm <- opop$mfmfmfmm[match(opop$mom, opop$pid)]
  opop$mfmfmffmm <- opop$mfmfmffm[match(opop$mom, opop$pid)]
  opop$mfmffmmmm <- opop$mfmffmmm[match(opop$mom, opop$pid)]
  opop$mfmffmfmm <- opop$mfmffmfm[match(opop$mom, opop$pid)]
  opop$mfmfffmmm <- opop$mfmfffmm[match(opop$mom, opop$pid)]
  opop$mfmffffmm <- opop$mfmffffm[match(opop$mom, opop$pid)]
  opop$mffmmmmmm <- opop$mffmmmmm[match(opop$mom, opop$pid)]
  opop$mffmmmfmm <- opop$mffmmmfm[match(opop$mom, opop$pid)]
  opop$mffmmfmmm <- opop$mffmmfmm[match(opop$mom, opop$pid)]
  opop$mffmmffmm <- opop$mffmmffm[match(opop$mom, opop$pid)]
  opop$mffmfmmmm <- opop$mffmfmmm[match(opop$mom, opop$pid)]
  opop$mffmfmfmm <- opop$mffmfmfm[match(opop$mom, opop$pid)]
  opop$mffmffmmm <- opop$mffmffmm[match(opop$mom, opop$pid)]
  opop$mffmfffmm <- opop$mffmfffm[match(opop$mom, opop$pid)]
  opop$mfffmmmmm <- opop$mfffmmmm[match(opop$mom, opop$pid)]
  opop$mfffmmfmm <- opop$mfffmmfm[match(opop$mom, opop$pid)]
  opop$mfffmfmmm <- opop$mfffmfmm[match(opop$mom, opop$pid)]
  opop$mfffmffmm <- opop$mfffmffm[match(opop$mom, opop$pid)]
  opop$mffffmmmm <- opop$mffffmmm[match(opop$mom, opop$pid)]
  opop$mffffmfmm <- opop$mffffmfm[match(opop$mom, opop$pid)]
  opop$mfffffmmm <- opop$mfffffmm[match(opop$mom, opop$pid)]
  opop$mffffffmm <- opop$mffffffm[match(opop$mom, opop$pid)]
  opop$mmmmmmmfm <- opop$mmmmmmmf[match(opop$mom, opop$pid)]
  opop$mmmmmmffm <- opop$mmmmmmff[match(opop$mom, opop$pid)]
  opop$mmmmmfmfm <- opop$mmmmmfmf[match(opop$mom, opop$pid)]
  opop$mmmmmfffm <- opop$mmmmmfff[match(opop$mom, opop$pid)]
  opop$mmmmfmmfm <- opop$mmmmfmmf[match(opop$mom, opop$pid)]
  opop$mmmmfmffm <- opop$mmmmfmff[match(opop$mom, opop$pid)]
  opop$mmmmffmfm <- opop$mmmmffmf[match(opop$mom, opop$pid)]
  opop$mmmmffffm <- opop$mmmmffff[match(opop$mom, opop$pid)]
  opop$mmmfmmmfm <- opop$mmmfmmmf[match(opop$mom, opop$pid)]
  opop$mmmfmmffm <- opop$mmmfmmff[match(opop$mom, opop$pid)]
  opop$mmmfmfmfm <- opop$mmmfmfmf[match(opop$mom, opop$pid)]
  opop$mmmfmfffm <- opop$mmmfmfff[match(opop$mom, opop$pid)]
  opop$mmmffmmfm <- opop$mmmffmmf[match(opop$mom, opop$pid)]
  opop$mmmffmffm <- opop$mmmffmff[match(opop$mom, opop$pid)]
  opop$mmmfffmfm <- opop$mmmfffmf[match(opop$mom, opop$pid)]
  opop$mmmfffffm <- opop$mmmfffff[match(opop$mom, opop$pid)]
  opop$mmfmmmmfm <- opop$mmfmmmmf[match(opop$mom, opop$pid)]
  opop$mmfmmmffm <- opop$mmfmmmff[match(opop$mom, opop$pid)]
  opop$mmfmmfmfm <- opop$mmfmmfmf[match(opop$mom, opop$pid)]
  opop$mmfmmfffm <- opop$mmfmmfff[match(opop$mom, opop$pid)]
  opop$mmfmfmmfm <- opop$mmfmfmmf[match(opop$mom, opop$pid)]
  opop$mmfmfmffm <- opop$mmfmfmff[match(opop$mom, opop$pid)]
  opop$mmfmffmfm <- opop$mmfmffmf[match(opop$mom, opop$pid)]
  opop$mmfmffffm <- opop$mmfmffff[match(opop$mom, opop$pid)]
  opop$mmffmmmfm <- opop$mmffmmmf[match(opop$mom, opop$pid)]
  opop$mmffmmffm <- opop$mmffmmff[match(opop$mom, opop$pid)]
  opop$mmffmfmfm <- opop$mmffmfmf[match(opop$mom, opop$pid)]
  opop$mmffmfffm <- opop$mmffmfff[match(opop$mom, opop$pid)]
  opop$mmfffmmfm <- opop$mmfffmmf[match(opop$mom, opop$pid)]
  opop$mmfffmffm <- opop$mmfffmff[match(opop$mom, opop$pid)]
  opop$mmffffmfm <- opop$mmffffmf[match(opop$mom, opop$pid)]
  opop$mmffffffm <- opop$mmffffff[match(opop$mom, opop$pid)]
  opop$mfmmmmmfm <- opop$mfmmmmmf[match(opop$mom, opop$pid)]
  opop$mfmmmmffm <- opop$mfmmmmff[match(opop$mom, opop$pid)]
  opop$mfmmmfmfm <- opop$mfmmmfmf[match(opop$mom, opop$pid)]
  opop$mfmmmfffm <- opop$mfmmmfff[match(opop$mom, opop$pid)]
  opop$mfmmfmmfm <- opop$mfmmfmmf[match(opop$mom, opop$pid)]
  opop$mfmmfmffm <- opop$mfmmfmff[match(opop$mom, opop$pid)]
  opop$mfmmffmfm <- opop$mfmmffmf[match(opop$mom, opop$pid)]
  opop$mfmmffffm <- opop$mfmmffff[match(opop$mom, opop$pid)]
  opop$mfmfmmmfm <- opop$mfmfmmmf[match(opop$mom, opop$pid)]
  opop$mfmfmmffm <- opop$mfmfmmff[match(opop$mom, opop$pid)]
  opop$mfmfmfmfm <- opop$mfmfmfmf[match(opop$mom, opop$pid)]
  opop$mfmfmfffm <- opop$mfmfmfff[match(opop$mom, opop$pid)]
  opop$mfmffmmfm <- opop$mfmffmmf[match(opop$mom, opop$pid)]
  opop$mfmffmffm <- opop$mfmffmff[match(opop$mom, opop$pid)]
  opop$mfmfffmfm <- opop$mfmfffmf[match(opop$mom, opop$pid)]
  opop$mfmfffffm <- opop$mfmfffff[match(opop$mom, opop$pid)]
  opop$mffmmmmfm <- opop$mffmmmmf[match(opop$mom, opop$pid)]
  opop$mffmmmffm <- opop$mffmmmff[match(opop$mom, opop$pid)]
  opop$mffmmfmfm <- opop$mffmmfmf[match(opop$mom, opop$pid)]
  opop$mffmmfffm <- opop$mffmmfff[match(opop$mom, opop$pid)]
  opop$mffmfmmfm <- opop$mffmfmmf[match(opop$mom, opop$pid)]
  opop$mffmfmffm <- opop$mffmfmff[match(opop$mom, opop$pid)]
  opop$mffmffmfm <- opop$mffmffmf[match(opop$mom, opop$pid)]
  opop$mffmffffm <- opop$mffmffff[match(opop$mom, opop$pid)]
  opop$mfffmmmfm <- opop$mfffmmmf[match(opop$mom, opop$pid)]
  opop$mfffmmffm <- opop$mfffmmff[match(opop$mom, opop$pid)]
  opop$mfffmfmfm <- opop$mfffmfmf[match(opop$mom, opop$pid)]
  opop$mfffmfffm <- opop$mfffmfff[match(opop$mom, opop$pid)]
  opop$mffffmmfm <- opop$mffffmmf[match(opop$mom, opop$pid)]
  opop$mffffmffm <- opop$mffffmff[match(opop$mom, opop$pid)]
  opop$mfffffmfm <- opop$mfffffmf[match(opop$mom, opop$pid)]
  opop$mfffffffm <- opop$mfffffff[match(opop$mom, opop$pid)]
  opop$mmmmmmmmf <- opop$mmmmmmmm[match(opop$pop, opop$pid)]
  opop$mmmmmmfmf <- opop$mmmmmmfm[match(opop$pop, opop$pid)]
  opop$mmmmmfmmf <- opop$mmmmmfmm[match(opop$pop, opop$pid)]
  opop$mmmmmffmf <- opop$mmmmmffm[match(opop$pop, opop$pid)]
  opop$mmmmfmmmf <- opop$mmmmfmmm[match(opop$pop, opop$pid)]
  opop$mmmmfmfmf <- opop$mmmmfmfm[match(opop$pop, opop$pid)]
  opop$mmmmffmmf <- opop$mmmmffmm[match(opop$pop, opop$pid)]
  opop$mmmmfffmf <- opop$mmmmfffm[match(opop$pop, opop$pid)]
  opop$mmmfmmmmf <- opop$mmmfmmmm[match(opop$pop, opop$pid)]
  opop$mmmfmmfmf <- opop$mmmfmmfm[match(opop$pop, opop$pid)]
  opop$mmmfmfmmf <- opop$mmmfmfmm[match(opop$pop, opop$pid)]
  opop$mmmfmffmf <- opop$mmmfmffm[match(opop$pop, opop$pid)]
  opop$mmmffmmmf <- opop$mmmffmmm[match(opop$pop, opop$pid)]
  opop$mmmffmfmf <- opop$mmmffmfm[match(opop$pop, opop$pid)]
  opop$mmmfffmmf <- opop$mmmfffmm[match(opop$pop, opop$pid)]
  opop$mmmffffmf <- opop$mmmffffm[match(opop$pop, opop$pid)]
  opop$mmfmmmmmf <- opop$mmfmmmmm[match(opop$pop, opop$pid)]
  opop$mmfmmmfmf <- opop$mmfmmmfm[match(opop$pop, opop$pid)]
  opop$mmfmmfmmf <- opop$mmfmmfmm[match(opop$pop, opop$pid)]
  opop$mmfmmffmf <- opop$mmfmmffm[match(opop$pop, opop$pid)]
  opop$mmfmfmmmf <- opop$mmfmfmmm[match(opop$pop, opop$pid)]
  opop$mmfmfmfmf <- opop$mmfmfmfm[match(opop$pop, opop$pid)]
  opop$mmfmffmmf <- opop$mmfmffmm[match(opop$pop, opop$pid)]
  opop$mmfmfffmf <- opop$mmfmfffm[match(opop$pop, opop$pid)]
  opop$mmffmmmmf <- opop$mmffmmmm[match(opop$pop, opop$pid)]
  opop$mmffmmfmf <- opop$mmffmmfm[match(opop$pop, opop$pid)]
  opop$mmffmfmmf <- opop$mmffmfmm[match(opop$pop, opop$pid)]
  opop$mmffmffmf <- opop$mmffmffm[match(opop$pop, opop$pid)]
  opop$mmfffmmmf <- opop$mmfffmmm[match(opop$pop, opop$pid)]
  opop$mmfffmfmf <- opop$mmfffmfm[match(opop$pop, opop$pid)]
  opop$mmffffmmf <- opop$mmffffmm[match(opop$pop, opop$pid)]
  opop$mmfffffmf <- opop$mmfffffm[match(opop$pop, opop$pid)]
  opop$mfmmmmmmf <- opop$mfmmmmmm[match(opop$pop, opop$pid)]
  opop$mfmmmmfmf <- opop$mfmmmmfm[match(opop$pop, opop$pid)]
  opop$mfmmmfmmf <- opop$mfmmmfmm[match(opop$pop, opop$pid)]
  opop$mfmmmffmf <- opop$mfmmmffm[match(opop$pop, opop$pid)]
  opop$mfmmfmmmf <- opop$mfmmfmmm[match(opop$pop, opop$pid)]
  opop$mfmmfmfmf <- opop$mfmmfmfm[match(opop$pop, opop$pid)]
  opop$mfmmffmmf <- opop$mfmmffmm[match(opop$pop, opop$pid)]
  opop$mfmmfffmf <- opop$mfmmfffm[match(opop$pop, opop$pid)]
  opop$mfmfmmmmf <- opop$mfmfmmmm[match(opop$pop, opop$pid)]
  opop$mfmfmmfmf <- opop$mfmfmmfm[match(opop$pop, opop$pid)]
  opop$mfmfmfmmf <- opop$mfmfmfmm[match(opop$pop, opop$pid)]
  opop$mfmfmffmf <- opop$mfmfmffm[match(opop$pop, opop$pid)]
  opop$mfmffmmmf <- opop$mfmffmmm[match(opop$pop, opop$pid)]
  opop$mfmffmfmf <- opop$mfmffmfm[match(opop$pop, opop$pid)]
  opop$mfmfffmmf <- opop$mfmfffmm[match(opop$pop, opop$pid)]
  opop$mfmffffmf <- opop$mfmffffm[match(opop$pop, opop$pid)]
  opop$mffmmmmmf <- opop$mffmmmmm[match(opop$pop, opop$pid)]
  opop$mffmmmfmf <- opop$mffmmmfm[match(opop$pop, opop$pid)]
  opop$mffmmfmmf <- opop$mffmmfmm[match(opop$pop, opop$pid)]
  opop$mffmmffmf <- opop$mffmmffm[match(opop$pop, opop$pid)]
  opop$mffmfmmmf <- opop$mffmfmmm[match(opop$pop, opop$pid)]
  opop$mffmfmfmf <- opop$mffmfmfm[match(opop$pop, opop$pid)]
  opop$mffmffmmf <- opop$mffmffmm[match(opop$pop, opop$pid)]
  opop$mffmfffmf <- opop$mffmfffm[match(opop$pop, opop$pid)]
  opop$mfffmmmmf <- opop$mfffmmmm[match(opop$pop, opop$pid)]
  opop$mfffmmfmf <- opop$mfffmmfm[match(opop$pop, opop$pid)]
  opop$mfffmfmmf <- opop$mfffmfmm[match(opop$pop, opop$pid)]
  opop$mfffmffmf <- opop$mfffmffm[match(opop$pop, opop$pid)]
  opop$mffffmmmf <- opop$mffffmmm[match(opop$pop, opop$pid)]
  opop$mffffmfmf <- opop$mffffmfm[match(opop$pop, opop$pid)]
  opop$mfffffmmf <- opop$mfffffmm[match(opop$pop, opop$pid)]
  opop$mffffffmf <- opop$mffffffm[match(opop$pop, opop$pid)]
  opop$mmmmmmmff <- opop$mmmmmmmf[match(opop$pop, opop$pid)]
  opop$mmmmmmfff <- opop$mmmmmmff[match(opop$pop, opop$pid)]
  opop$mmmmmfmff <- opop$mmmmmfmf[match(opop$pop, opop$pid)]
  opop$mmmmmffff <- opop$mmmmmfff[match(opop$pop, opop$pid)]
  opop$mmmmfmmff <- opop$mmmmfmmf[match(opop$pop, opop$pid)]
  opop$mmmmfmfff <- opop$mmmmfmff[match(opop$pop, opop$pid)]
  opop$mmmmffmff <- opop$mmmmffmf[match(opop$pop, opop$pid)]
  opop$mmmmfffff <- opop$mmmmffff[match(opop$pop, opop$pid)]
  opop$mmmfmmmff <- opop$mmmfmmmf[match(opop$pop, opop$pid)]
  opop$mmmfmmfff <- opop$mmmfmmff[match(opop$pop, opop$pid)]
  opop$mmmfmfmff <- opop$mmmfmfmf[match(opop$pop, opop$pid)]
  opop$mmmfmffff <- opop$mmmfmfff[match(opop$pop, opop$pid)]
  opop$mmmffmmff <- opop$mmmffmmf[match(opop$pop, opop$pid)]
  opop$mmmffmfff <- opop$mmmffmff[match(opop$pop, opop$pid)]
  opop$mmmfffmff <- opop$mmmfffmf[match(opop$pop, opop$pid)]
  opop$mmmffffff <- opop$mmmfffff[match(opop$pop, opop$pid)]
  opop$mmfmmmmff <- opop$mmfmmmmf[match(opop$pop, opop$pid)]
  opop$mmfmmmfff <- opop$mmfmmmff[match(opop$pop, opop$pid)]
  opop$mmfmmfmff <- opop$mmfmmfmf[match(opop$pop, opop$pid)]
  opop$mmfmmffff <- opop$mmfmmfff[match(opop$pop, opop$pid)]
  opop$mmfmfmmff <- opop$mmfmfmmf[match(opop$pop, opop$pid)]
  opop$mmfmfmfff <- opop$mmfmfmff[match(opop$pop, opop$pid)]
  opop$mmfmffmff <- opop$mmfmffmf[match(opop$pop, opop$pid)]
  opop$mmfmfffff <- opop$mmfmffff[match(opop$pop, opop$pid)]
  opop$mmffmmmff <- opop$mmffmmmf[match(opop$pop, opop$pid)]
  opop$mmffmmfff <- opop$mmffmmff[match(opop$pop, opop$pid)]
  opop$mmffmfmff <- opop$mmffmfmf[match(opop$pop, opop$pid)]
  opop$mmffmffff <- opop$mmffmfff[match(opop$pop, opop$pid)]
  opop$mmfffmmff <- opop$mmfffmmf[match(opop$pop, opop$pid)]
  opop$mmfffmfff <- opop$mmfffmff[match(opop$pop, opop$pid)]
  opop$mmffffmff <- opop$mmffffmf[match(opop$pop, opop$pid)]
  opop$mmfffffff <- opop$mmffffff[match(opop$pop, opop$pid)]
  opop$mfmmmmmff <- opop$mfmmmmmf[match(opop$pop, opop$pid)]
  opop$mfmmmmfff <- opop$mfmmmmff[match(opop$pop, opop$pid)]
  opop$mfmmmfmff <- opop$mfmmmfmf[match(opop$pop, opop$pid)]
  opop$mfmmmffff <- opop$mfmmmfff[match(opop$pop, opop$pid)]
  opop$mfmmfmmff <- opop$mfmmfmmf[match(opop$pop, opop$pid)]
  opop$mfmmfmfff <- opop$mfmmfmff[match(opop$pop, opop$pid)]
  opop$mfmmffmff <- opop$mfmmffmf[match(opop$pop, opop$pid)]
  opop$mfmmfffff <- opop$mfmmffff[match(opop$pop, opop$pid)]
  opop$mfmfmmmff <- opop$mfmfmmmf[match(opop$pop, opop$pid)]
  opop$mfmfmmfff <- opop$mfmfmmff[match(opop$pop, opop$pid)]
  opop$mfmfmfmff <- opop$mfmfmfmf[match(opop$pop, opop$pid)]
  opop$mfmfmffff <- opop$mfmfmfff[match(opop$pop, opop$pid)]
  opop$mfmffmmff <- opop$mfmffmmf[match(opop$pop, opop$pid)]
  opop$mfmffmfff <- opop$mfmffmff[match(opop$pop, opop$pid)]
  opop$mfmfffmff <- opop$mfmfffmf[match(opop$pop, opop$pid)]
  opop$mfmffffff <- opop$mfmfffff[match(opop$pop, opop$pid)]
  opop$mffmmmmff <- opop$mffmmmmf[match(opop$pop, opop$pid)]
  opop$mffmmmfff <- opop$mffmmmff[match(opop$pop, opop$pid)]
  opop$mffmmfmff <- opop$mffmmfmf[match(opop$pop, opop$pid)]
  opop$mffmmffff <- opop$mffmmfff[match(opop$pop, opop$pid)]
  opop$mffmfmmff <- opop$mffmfmmf[match(opop$pop, opop$pid)]
  opop$mffmfmfff <- opop$mffmfmff[match(opop$pop, opop$pid)]
  opop$mffmffmff <- opop$mffmffmf[match(opop$pop, opop$pid)]
  opop$mffmfffff <- opop$mffmffff[match(opop$pop, opop$pid)]
  opop$mfffmmmff <- opop$mfffmmmf[match(opop$pop, opop$pid)]
  opop$mfffmmfff <- opop$mfffmmff[match(opop$pop, opop$pid)]
  opop$mfffmfmff <- opop$mfffmfmf[match(opop$pop, opop$pid)]
  opop$mfffmffff <- opop$mfffmfff[match(opop$pop, opop$pid)]
  opop$mffffmmff <- opop$mffffmmf[match(opop$pop, opop$pid)]
  opop$mffffmfff <- opop$mffffmff[match(opop$pop, opop$pid)]
  opop$mfffffmff <- opop$mfffffmf[match(opop$pop, opop$pid)]
  opop$mffffffff <- opop$mfffffff[match(opop$pop, opop$pid)]
  opop$fmmmmmmmm <- opop$fmmmmmmm[match(opop$mom, opop$pid)]
  opop$fmmmmmfmm <- opop$fmmmmmfm[match(opop$mom, opop$pid)]
  opop$fmmmmfmmm <- opop$fmmmmfmm[match(opop$mom, opop$pid)]
  opop$fmmmmffmm <- opop$fmmmmffm[match(opop$mom, opop$pid)]
  opop$fmmmfmmmm <- opop$fmmmfmmm[match(opop$mom, opop$pid)]
  opop$fmmmfmfmm <- opop$fmmmfmfm[match(opop$mom, opop$pid)]
  opop$fmmmffmmm <- opop$fmmmffmm[match(opop$mom, opop$pid)]
  opop$fmmmfffmm <- opop$fmmmfffm[match(opop$mom, opop$pid)]
  opop$fmmfmmmmm <- opop$fmmfmmmm[match(opop$mom, opop$pid)]
  opop$fmmfmmfmm <- opop$fmmfmmfm[match(opop$mom, opop$pid)]
  opop$fmmfmfmmm <- opop$fmmfmfmm[match(opop$mom, opop$pid)]
  opop$fmmfmffmm <- opop$fmmfmffm[match(opop$mom, opop$pid)]
  opop$fmmffmmmm <- opop$fmmffmmm[match(opop$mom, opop$pid)]
  opop$fmmffmfmm <- opop$fmmffmfm[match(opop$mom, opop$pid)]
  opop$fmmfffmmm <- opop$fmmfffmm[match(opop$mom, opop$pid)]
  opop$fmmffffmm <- opop$fmmffffm[match(opop$mom, opop$pid)]
  opop$fmfmmmmmm <- opop$fmfmmmmm[match(opop$mom, opop$pid)]
  opop$fmfmmmfmm <- opop$fmfmmmfm[match(opop$mom, opop$pid)]
  opop$fmfmmfmmm <- opop$fmfmmfmm[match(opop$mom, opop$pid)]
  opop$fmfmmffmm <- opop$fmfmmffm[match(opop$mom, opop$pid)]
  opop$fmfmfmmmm <- opop$fmfmfmmm[match(opop$mom, opop$pid)]
  opop$fmfmfmfmm <- opop$fmfmfmfm[match(opop$mom, opop$pid)]
  opop$fmfmffmmm <- opop$fmfmffmm[match(opop$mom, opop$pid)]
  opop$fmfmfffmm <- opop$fmfmfffm[match(opop$mom, opop$pid)]
  opop$fmffmmmmm <- opop$fmffmmmm[match(opop$mom, opop$pid)]
  opop$fmffmmfmm <- opop$fmffmmfm[match(opop$mom, opop$pid)]
  opop$fmffmfmmm <- opop$fmffmfmm[match(opop$mom, opop$pid)]
  opop$fmffmffmm <- opop$fmffmffm[match(opop$mom, opop$pid)]
  opop$fmfffmmmm <- opop$fmfffmmm[match(opop$mom, opop$pid)]
  opop$fmfffmfmm <- opop$fmfffmfm[match(opop$mom, opop$pid)]
  opop$fmffffmmm <- opop$fmffffmm[match(opop$mom, opop$pid)]
  opop$fmfffffmm <- opop$fmfffffm[match(opop$mom, opop$pid)]
  opop$ffmmmmmmm <- opop$ffmmmmmm[match(opop$mom, opop$pid)]
  opop$ffmmmmfmm <- opop$ffmmmmfm[match(opop$mom, opop$pid)]
  opop$ffmmmfmmm <- opop$ffmmmfmm[match(opop$mom, opop$pid)]
  opop$ffmmmffmm <- opop$ffmmmffm[match(opop$mom, opop$pid)]
  opop$ffmmfmmmm <- opop$ffmmfmmm[match(opop$mom, opop$pid)]
  opop$ffmmfmfmm <- opop$ffmmfmfm[match(opop$mom, opop$pid)]
  opop$ffmmffmmm <- opop$ffmmffmm[match(opop$mom, opop$pid)]
  opop$ffmmfffmm <- opop$ffmmfffm[match(opop$mom, opop$pid)]
  opop$ffmfmmmmm <- opop$ffmfmmmm[match(opop$mom, opop$pid)]
  opop$ffmfmmfmm <- opop$ffmfmmfm[match(opop$mom, opop$pid)]
  opop$ffmfmfmmm <- opop$ffmfmfmm[match(opop$mom, opop$pid)]
  opop$ffmfmffmm <- opop$ffmfmffm[match(opop$mom, opop$pid)]
  opop$ffmffmmmm <- opop$ffmffmmm[match(opop$mom, opop$pid)]
  opop$ffmffmfmm <- opop$ffmffmfm[match(opop$mom, opop$pid)]
  opop$ffmfffmmm <- opop$ffmfffmm[match(opop$mom, opop$pid)]
  opop$ffmffffmm <- opop$ffmffffm[match(opop$mom, opop$pid)]
  opop$fffmmmmmm <- opop$fffmmmmm[match(opop$mom, opop$pid)]
  opop$fffmmmfmm <- opop$fffmmmfm[match(opop$mom, opop$pid)]
  opop$fffmmfmmm <- opop$fffmmfmm[match(opop$mom, opop$pid)]
  opop$fffmmffmm <- opop$fffmmffm[match(opop$mom, opop$pid)]
  opop$fffmfmmmm <- opop$fffmfmmm[match(opop$mom, opop$pid)]
  opop$fffmfmfmm <- opop$fffmfmfm[match(opop$mom, opop$pid)]
  opop$fffmffmmm <- opop$fffmffmm[match(opop$mom, opop$pid)]
  opop$fffmfffmm <- opop$fffmfffm[match(opop$mom, opop$pid)]
  opop$ffffmmmmm <- opop$ffffmmmm[match(opop$mom, opop$pid)]
  opop$ffffmmfmm <- opop$ffffmmfm[match(opop$mom, opop$pid)]
  opop$ffffmfmmm <- opop$ffffmfmm[match(opop$mom, opop$pid)]
  opop$ffffmffmm <- opop$ffffmffm[match(opop$mom, opop$pid)]
  opop$fffffmmmm <- opop$fffffmmm[match(opop$mom, opop$pid)]
  opop$fffffmfmm <- opop$fffffmfm[match(opop$mom, opop$pid)]
  opop$ffffffmmm <- opop$ffffffmm[match(opop$mom, opop$pid)]
  opop$fffffffmm <- opop$fffffffm[match(opop$mom, opop$pid)]
  opop$fmmmmmmfm <- opop$fmmmmmmf[match(opop$mom, opop$pid)]
  opop$fmmmmmffm <- opop$fmmmmmff[match(opop$mom, opop$pid)]
  opop$fmmmmfmfm <- opop$fmmmmfmf[match(opop$mom, opop$pid)]
  opop$fmmmmfffm <- opop$fmmmmfff[match(opop$mom, opop$pid)]
  opop$fmmmfmmfm <- opop$fmmmfmmf[match(opop$mom, opop$pid)]
  opop$fmmmfmffm <- opop$fmmmfmff[match(opop$mom, opop$pid)]
  opop$fmmmffmfm <- opop$fmmmffmf[match(opop$mom, opop$pid)]
  opop$fmmmffffm <- opop$fmmmffff[match(opop$mom, opop$pid)]
  opop$fmmfmmmfm <- opop$fmmfmmmf[match(opop$mom, opop$pid)]
  opop$fmmfmmffm <- opop$fmmfmmff[match(opop$mom, opop$pid)]
  opop$fmmfmfmfm <- opop$fmmfmfmf[match(opop$mom, opop$pid)]
  opop$fmmfmfffm <- opop$fmmfmfff[match(opop$mom, opop$pid)]
  opop$fmmffmmfm <- opop$fmmffmmf[match(opop$mom, opop$pid)]
  opop$fmmffmffm <- opop$fmmffmff[match(opop$mom, opop$pid)]
  opop$fmmfffmfm <- opop$fmmfffmf[match(opop$mom, opop$pid)]
  opop$fmmfffffm <- opop$fmmfffff[match(opop$mom, opop$pid)]
  opop$fmfmmmmfm <- opop$fmfmmmmf[match(opop$mom, opop$pid)]
  opop$fmfmmmffm <- opop$fmfmmmff[match(opop$mom, opop$pid)]
  opop$fmfmmfmfm <- opop$fmfmmfmf[match(opop$mom, opop$pid)]
  opop$fmfmmfffm <- opop$fmfmmfff[match(opop$mom, opop$pid)]
  opop$fmfmfmmfm <- opop$fmfmfmmf[match(opop$mom, opop$pid)]
  opop$fmfmfmffm <- opop$fmfmfmff[match(opop$mom, opop$pid)]
  opop$fmfmffmfm <- opop$fmfmffmf[match(opop$mom, opop$pid)]
  opop$fmfmffffm <- opop$fmfmffff[match(opop$mom, opop$pid)]
  opop$fmffmmmfm <- opop$fmffmmmf[match(opop$mom, opop$pid)]
  opop$fmffmmffm <- opop$fmffmmff[match(opop$mom, opop$pid)]
  opop$fmffmfmfm <- opop$fmffmfmf[match(opop$mom, opop$pid)]
  opop$fmffmfffm <- opop$fmffmfff[match(opop$mom, opop$pid)]
  opop$fmfffmmfm <- opop$fmfffmmf[match(opop$mom, opop$pid)]
  opop$fmfffmffm <- opop$fmfffmff[match(opop$mom, opop$pid)]
  opop$fmffffmfm <- opop$fmffffmf[match(opop$mom, opop$pid)]
  opop$fmffffffm <- opop$fmffffff[match(opop$mom, opop$pid)]
  opop$ffmmmmmfm <- opop$ffmmmmmf[match(opop$mom, opop$pid)]
  opop$ffmmmmffm <- opop$ffmmmmff[match(opop$mom, opop$pid)]
  opop$ffmmmfmfm <- opop$ffmmmfmf[match(opop$mom, opop$pid)]
  opop$ffmmmfffm <- opop$ffmmmfff[match(opop$mom, opop$pid)]
  opop$ffmmfmmfm <- opop$ffmmfmmf[match(opop$mom, opop$pid)]
  opop$ffmmfmffm <- opop$ffmmfmff[match(opop$mom, opop$pid)]
  opop$ffmmffmfm <- opop$ffmmffmf[match(opop$mom, opop$pid)]
  opop$ffmmffffm <- opop$ffmmffff[match(opop$mom, opop$pid)]
  opop$ffmfmmmfm <- opop$ffmfmmmf[match(opop$mom, opop$pid)]
  opop$ffmfmmffm <- opop$ffmfmmff[match(opop$mom, opop$pid)]
  opop$ffmfmfmfm <- opop$ffmfmfmf[match(opop$mom, opop$pid)]
  opop$ffmfmfffm <- opop$ffmfmfff[match(opop$mom, opop$pid)]
  opop$ffmffmmfm <- opop$ffmffmmf[match(opop$mom, opop$pid)]
  opop$ffmffmffm <- opop$ffmffmff[match(opop$mom, opop$pid)]
  opop$ffmfffmfm <- opop$ffmfffmf[match(opop$mom, opop$pid)]
  opop$ffmfffffm <- opop$ffmfffff[match(opop$mom, opop$pid)]
  opop$fffmmmmfm <- opop$fffmmmmf[match(opop$mom, opop$pid)]
  opop$fffmmmffm <- opop$fffmmmff[match(opop$mom, opop$pid)]
  opop$fffmmfmfm <- opop$fffmmfmf[match(opop$mom, opop$pid)]
  opop$fffmmfffm <- opop$fffmmfff[match(opop$mom, opop$pid)]
  opop$fffmfmmfm <- opop$fffmfmmf[match(opop$mom, opop$pid)]
  opop$fffmfmffm <- opop$fffmfmff[match(opop$mom, opop$pid)]
  opop$fffmffmfm <- opop$fffmffmf[match(opop$mom, opop$pid)]
  opop$fffmffffm <- opop$fffmffff[match(opop$mom, opop$pid)]
  opop$ffffmmmfm <- opop$ffffmmmf[match(opop$mom, opop$pid)]
  opop$ffffmmffm <- opop$ffffmmff[match(opop$mom, opop$pid)]
  opop$ffffmfmfm <- opop$ffffmfmf[match(opop$mom, opop$pid)]
  opop$ffffmfffm <- opop$ffffmfff[match(opop$mom, opop$pid)]
  opop$fffffmmfm <- opop$fffffmmf[match(opop$mom, opop$pid)]
  opop$fffffmffm <- opop$fffffmff[match(opop$mom, opop$pid)]
  opop$ffffffmfm <- opop$ffffffmf[match(opop$mom, opop$pid)]
  opop$ffffffffm <- opop$ffffffff[match(opop$mom, opop$pid)]
  opop$fmmmmmmmf <- opop$fmmmmmmm[match(opop$pop, opop$pid)]
  opop$fmmmmmfmf <- opop$fmmmmmfm[match(opop$pop, opop$pid)]
  opop$fmmmmfmmf <- opop$fmmmmfmm[match(opop$pop, opop$pid)]
  opop$fmmmmffmf <- opop$fmmmmffm[match(opop$pop, opop$pid)]
  opop$fmmmfmmmf <- opop$fmmmfmmm[match(opop$pop, opop$pid)]
  opop$fmmmfmfmf <- opop$fmmmfmfm[match(opop$pop, opop$pid)]
  opop$fmmmffmmf <- opop$fmmmffmm[match(opop$pop, opop$pid)]
  opop$fmmmfffmf <- opop$fmmmfffm[match(opop$pop, opop$pid)]
  opop$fmmfmmmmf <- opop$fmmfmmmm[match(opop$pop, opop$pid)]
  opop$fmmfmmfmf <- opop$fmmfmmfm[match(opop$pop, opop$pid)]
  opop$fmmfmfmmf <- opop$fmmfmfmm[match(opop$pop, opop$pid)]
  opop$fmmfmffmf <- opop$fmmfmffm[match(opop$pop, opop$pid)]
  opop$fmmffmmmf <- opop$fmmffmmm[match(opop$pop, opop$pid)]
  opop$fmmffmfmf <- opop$fmmffmfm[match(opop$pop, opop$pid)]
  opop$fmmfffmmf <- opop$fmmfffmm[match(opop$pop, opop$pid)]
  opop$fmmffffmf <- opop$fmmffffm[match(opop$pop, opop$pid)]
  opop$fmfmmmmmf <- opop$fmfmmmmm[match(opop$pop, opop$pid)]
  opop$fmfmmmfmf <- opop$fmfmmmfm[match(opop$pop, opop$pid)]
  opop$fmfmmfmmf <- opop$fmfmmfmm[match(opop$pop, opop$pid)]
  opop$fmfmmffmf <- opop$fmfmmffm[match(opop$pop, opop$pid)]
  opop$fmfmfmmmf <- opop$fmfmfmmm[match(opop$pop, opop$pid)]
  opop$fmfmfmfmf <- opop$fmfmfmfm[match(opop$pop, opop$pid)]
  opop$fmfmffmmf <- opop$fmfmffmm[match(opop$pop, opop$pid)]
  opop$fmfmfffmf <- opop$fmfmfffm[match(opop$pop, opop$pid)]
  opop$fmffmmmmf <- opop$fmffmmmm[match(opop$pop, opop$pid)]
  opop$fmffmmfmf <- opop$fmffmmfm[match(opop$pop, opop$pid)]
  opop$fmffmfmmf <- opop$fmffmfmm[match(opop$pop, opop$pid)]
  opop$fmffmffmf <- opop$fmffmffm[match(opop$pop, opop$pid)]
  opop$fmfffmmmf <- opop$fmfffmmm[match(opop$pop, opop$pid)]
  opop$fmfffmfmf <- opop$fmfffmfm[match(opop$pop, opop$pid)]
  opop$fmffffmmf <- opop$fmffffmm[match(opop$pop, opop$pid)]
  opop$fmfffffmf <- opop$fmfffffm[match(opop$pop, opop$pid)]
  opop$ffmmmmmmf <- opop$ffmmmmmm[match(opop$pop, opop$pid)]
  opop$ffmmmmfmf <- opop$ffmmmmfm[match(opop$pop, opop$pid)]
  opop$ffmmmfmmf <- opop$ffmmmfmm[match(opop$pop, opop$pid)]
  opop$ffmmmffmf <- opop$ffmmmffm[match(opop$pop, opop$pid)]
  opop$ffmmfmmmf <- opop$ffmmfmmm[match(opop$pop, opop$pid)]
  opop$ffmmfmfmf <- opop$ffmmfmfm[match(opop$pop, opop$pid)]
  opop$ffmmffmmf <- opop$ffmmffmm[match(opop$pop, opop$pid)]
  opop$ffmmfffmf <- opop$ffmmfffm[match(opop$pop, opop$pid)]
  opop$ffmfmmmmf <- opop$ffmfmmmm[match(opop$pop, opop$pid)]
  opop$ffmfmmfmf <- opop$ffmfmmfm[match(opop$pop, opop$pid)]
  opop$ffmfmfmmf <- opop$ffmfmfmm[match(opop$pop, opop$pid)]
  opop$ffmfmffmf <- opop$ffmfmffm[match(opop$pop, opop$pid)]
  opop$ffmffmmmf <- opop$ffmffmmm[match(opop$pop, opop$pid)]
  opop$ffmffmfmf <- opop$ffmffmfm[match(opop$pop, opop$pid)]
  opop$ffmfffmmf <- opop$ffmfffmm[match(opop$pop, opop$pid)]
  opop$ffmffffmf <- opop$ffmffffm[match(opop$pop, opop$pid)]
  opop$fffmmmmmf <- opop$fffmmmmm[match(opop$pop, opop$pid)]
  opop$fffmmmfmf <- opop$fffmmmfm[match(opop$pop, opop$pid)]
  opop$fffmmfmmf <- opop$fffmmfmm[match(opop$pop, opop$pid)]
  opop$fffmmffmf <- opop$fffmmffm[match(opop$pop, opop$pid)]
  opop$fffmfmmmf <- opop$fffmfmmm[match(opop$pop, opop$pid)]
  opop$fffmfmfmf <- opop$fffmfmfm[match(opop$pop, opop$pid)]
  opop$fffmffmmf <- opop$fffmffmm[match(opop$pop, opop$pid)]
  opop$fffmfffmf <- opop$fffmfffm[match(opop$pop, opop$pid)]
  opop$ffffmmmmf <- opop$ffffmmmm[match(opop$pop, opop$pid)]
  opop$ffffmmfmf <- opop$ffffmmfm[match(opop$pop, opop$pid)]
  opop$ffffmfmmf <- opop$ffffmfmm[match(opop$pop, opop$pid)]
  opop$ffffmffmf <- opop$ffffmffm[match(opop$pop, opop$pid)]
  opop$fffffmmmf <- opop$fffffmmm[match(opop$pop, opop$pid)]
  opop$fffffmfmf <- opop$fffffmfm[match(opop$pop, opop$pid)]
  opop$ffffffmmf <- opop$ffffffmm[match(opop$pop, opop$pid)]
  opop$fffffffmf <- opop$fffffffm[match(opop$pop, opop$pid)]
  opop$fmmmmmmff <- opop$fmmmmmmf[match(opop$pop, opop$pid)]
  opop$fmmmmmfff <- opop$fmmmmmff[match(opop$pop, opop$pid)]
  opop$fmmmmfmff <- opop$fmmmmfmf[match(opop$pop, opop$pid)]
  opop$fmmmmffff <- opop$fmmmmfff[match(opop$pop, opop$pid)]
  opop$fmmmfmmff <- opop$fmmmfmmf[match(opop$pop, opop$pid)]
  opop$fmmmfmfff <- opop$fmmmfmff[match(opop$pop, opop$pid)]
  opop$fmmmffmff <- opop$fmmmffmf[match(opop$pop, opop$pid)]
  opop$fmmmfffff <- opop$fmmmffff[match(opop$pop, opop$pid)]
  opop$fmmfmmmff <- opop$fmmfmmmf[match(opop$pop, opop$pid)]
  opop$fmmfmmfff <- opop$fmmfmmff[match(opop$pop, opop$pid)]
  opop$fmmfmfmff <- opop$fmmfmfmf[match(opop$pop, opop$pid)]
  opop$fmmfmffff <- opop$fmmfmfff[match(opop$pop, opop$pid)]
  opop$fmmffmmff <- opop$fmmffmmf[match(opop$pop, opop$pid)]
  opop$fmmffmfff <- opop$fmmffmff[match(opop$pop, opop$pid)]
  opop$fmmfffmff <- opop$fmmfffmf[match(opop$pop, opop$pid)]
  opop$fmmffffff <- opop$fmmfffff[match(opop$pop, opop$pid)]
  opop$fmfmmmmff <- opop$fmfmmmmf[match(opop$pop, opop$pid)]
  opop$fmfmmmfff <- opop$fmfmmmff[match(opop$pop, opop$pid)]
  opop$fmfmmfmff <- opop$fmfmmfmf[match(opop$pop, opop$pid)]
  opop$fmfmmffff <- opop$fmfmmfff[match(opop$pop, opop$pid)]
  opop$fmfmfmmff <- opop$fmfmfmmf[match(opop$pop, opop$pid)]
  opop$fmfmfmfff <- opop$fmfmfmff[match(opop$pop, opop$pid)]
  opop$fmfmffmff <- opop$fmfmffmf[match(opop$pop, opop$pid)]
  opop$fmfmfffff <- opop$fmfmffff[match(opop$pop, opop$pid)]
  opop$fmffmmmff <- opop$fmffmmmf[match(opop$pop, opop$pid)]
  opop$fmffmmfff <- opop$fmffmmff[match(opop$pop, opop$pid)]
  opop$fmffmfmff <- opop$fmffmfmf[match(opop$pop, opop$pid)]
  opop$fmffmffff <- opop$fmffmfff[match(opop$pop, opop$pid)]
  opop$fmfffmmff <- opop$fmfffmmf[match(opop$pop, opop$pid)]
  opop$fmfffmfff <- opop$fmfffmff[match(opop$pop, opop$pid)]
  opop$fmffffmff <- opop$fmffffmf[match(opop$pop, opop$pid)]
  opop$fmfffffff <- opop$fmffffff[match(opop$pop, opop$pid)]
  opop$ffmmmmmff <- opop$ffmmmmmf[match(opop$pop, opop$pid)]
  opop$ffmmmmfff <- opop$ffmmmmff[match(opop$pop, opop$pid)]
  opop$ffmmmfmff <- opop$ffmmmfmf[match(opop$pop, opop$pid)]
  opop$ffmmmffff <- opop$ffmmmfff[match(opop$pop, opop$pid)]
  opop$ffmmfmmff <- opop$ffmmfmmf[match(opop$pop, opop$pid)]
  opop$ffmmfmfff <- opop$ffmmfmff[match(opop$pop, opop$pid)]
  opop$ffmmffmff <- opop$ffmmffmf[match(opop$pop, opop$pid)]
  opop$ffmmfffff <- opop$ffmmffff[match(opop$pop, opop$pid)]
  opop$ffmfmmmff <- opop$ffmfmmmf[match(opop$pop, opop$pid)]
  opop$ffmfmmfff <- opop$ffmfmmff[match(opop$pop, opop$pid)]
  opop$ffmfmfmff <- opop$ffmfmfmf[match(opop$pop, opop$pid)]
  opop$ffmfmffff <- opop$ffmfmfff[match(opop$pop, opop$pid)]
  opop$ffmffmmff <- opop$ffmffmmf[match(opop$pop, opop$pid)]
  opop$ffmffmfff <- opop$ffmffmff[match(opop$pop, opop$pid)]
  opop$ffmfffmff <- opop$ffmfffmf[match(opop$pop, opop$pid)]
  opop$ffmffffff <- opop$ffmfffff[match(opop$pop, opop$pid)]
  opop$fffmmmmff <- opop$fffmmmmf[match(opop$pop, opop$pid)]
  opop$fffmmmfff <- opop$fffmmmff[match(opop$pop, opop$pid)]
  opop$fffmmfmff <- opop$fffmmfmf[match(opop$pop, opop$pid)]
  opop$fffmmffff <- opop$fffmmfff[match(opop$pop, opop$pid)]
  opop$fffmfmmff <- opop$fffmfmmf[match(opop$pop, opop$pid)]
  opop$fffmfmfff <- opop$fffmfmff[match(opop$pop, opop$pid)]
  opop$fffmffmff <- opop$fffmffmf[match(opop$pop, opop$pid)]
  opop$fffmfffff <- opop$fffmffff[match(opop$pop, opop$pid)]
  opop$ffffmmmff <- opop$ffffmmmf[match(opop$pop, opop$pid)]
  opop$ffffmmfff <- opop$ffffmmff[match(opop$pop, opop$pid)]
  opop$ffffmfmff <- opop$ffffmfmf[match(opop$pop, opop$pid)]
  opop$ffffmffff <- opop$ffffmfff[match(opop$pop, opop$pid)]
  opop$fffffmmff <- opop$fffffmmf[match(opop$pop, opop$pid)]
  opop$fffffmfff <- opop$fffffmff[match(opop$pop, opop$pid)]
  opop$ffffffmff <- opop$ffffffmf[match(opop$pop, opop$pid)]
  opop$fffffffff <- opop$ffffffff[match(opop$pop, opop$pid)]
  
  
  kidsOf <- with(opop, { c(tapply(pid, mom, c), tapply(pid, pop, c)) })
  kidsOf["0"] <- NULL
  kidsOf["0"] <- NULL
  KidsOf <- list()
  KidsOf[as.numeric(names(kidsOf))] <- kidsOf
  ko <- function(KidsOf = KidsOf, p) {
    lapply(p, function(x) {
      unique(as.vector(unlist(KidsOf[x])))
    })
  }
  
  so <- function(opop = opop, p) {
    lapply(p, function(p) {
      as.vector(unlist(opop[opop$pid %in% p, c("spouse")]))
    })
  }
  
  res <- list()
  
  res$parents <- lapply(pid, function(pid) {
    as.vector(unlist(opop[opop$pid == pid, c("mom", "pop")]))
  })
  
  res$gparents <- lapply(pid, function(pid) {
    as.vector(unlist(opop[opop$pid == pid, c("mm", "mf", "fm", "ff")]))
  })
  
  res$ggparents <- lapply(pid, function(pid) {
    as.vector(unlist(opop[opop$pid == pid, c("mmm", "mfm", "mmf", "mff", 
                                             "fmm", "ffm", "fmf", "fff")]))
  }) 
  
  res$gggparents <- lapply(pid, function(pid) {
    as.vector(unlist(opop[opop$pid == pid, c("mmmm", "mmfm", "mfmm", "mffm", "mmmf", "mmff", "mfmf", "mfff",
                                             "fmmm", "fmfm", "ffmm", "fffm", "fmmf", "fmff", "ffmf", "ffff")]))
  })
  
  res$ggggparents <- lapply(pid, function(pid) {
    as.vector(unlist(opop[opop$pid == pid, c("mmmmm", "mmmfm", "mmfmm", "mmffm", "mfmmm", "mfmfm", "mffmm","mfffm", 
                                             "mmmmf", "mmmff", "mmfmf", "mmfff", "mfmmf", "mfmff", "mffmf", "mffff",
                                             "fmmmm", "fmmfm", "fmfmm", "fmffm", "ffmmm", "ffmfm", "fffmm", "ffffm", 
                                             "fmmmf", "fmmff", "fmfmf", "fmfff", "ffmmf", "ffmff", "fffmf", "fffff")]))
  })
  
  res$gggggparents <- lapply(pid, function(pid) {
    as.vector(unlist(opop[opop$pid == pid, c("mmmmmm", "mmmmfm", "mmmfmm", "mmmffm", "mmfmmm", "mmfmfm", "mmffmm", "mmfffm",
                                             "mfmmmm", "mfmmfm", "mfmfmm", "mfmffm", "mffmmm", "mffmfm", "mfffmm", "mffffm",
                                             "mmmmmf", "mmmmff", "mmmfmf", "mmmfff", "mmfmmf", "mmfmff", "mmffmf", "mmffff",
                                             "mfmmmf", "mfmmff", "mfmfmf", "mfmfff", "mffmmf", "mffmff", "mfffmf", "mfffff",
                                             "fmmmmm", "fmmmfm", "fmmfmm", "fmmffm", "fmfmmm", "fmfmfm", "fmffmm", "fmfffm",
                                             "ffmmmm", "ffmmfm", "ffmfmm", "ffmffm", "fffmmm", "fffmfm", "ffffmm", "fffffm",
                                             "fmmmmf", "fmmmff", "fmmfmf", "fmmfff", "fmfmmf", "fmfmff", "fmffmf", "fmffff",
                                             "ffmmmf", "ffmmff", "ffmfmf", "ffmfff", "fffmmf", "fffmff", "ffffmf", "ffffff")]))
  })

  res$ggggggparents <- lapply(pid, function(pid) {
    as.vector(unlist(opop[opop$pid == pid, c("mmmmmmm", "mmmmmfm", "mmmmfmm", "mmmmffm", "mmmfmmm", "mmmfmfm", "mmmffmm", "mmmfffm",
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
                                             "fffmmmf", "fffmmff", "fffmfmf", "fffmfff", "ffffmmf", "ffffmff", "fffffmf", "fffffff")]))
  })
  
  res$gggggggparents <- lapply(pid, function(pid) {
    as.vector(unlist(opop[opop$pid == pid, c("mmmmmmmm", "mmmmmmfm", "mmmmmfmm", "mmmmmffm", "mmmmfmmm", "mmmmfmfm", "mmmmffmm", "mmmmfffm",
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
                                             "ffffmmmf", "ffffmmff", "ffffmfmf", "ffffmfff", "fffffmmf", "fffffmff", "ffffffmf", "ffffffff")]))
  })


  res$ggggggggparents <- lapply(pid, function(pid) {
    as.vector(unlist(opop[opop$pid == pid, c("mmmmmmmmm", "mmmmmmfmm", "mmmmmfmmm", "mmmmmffmm", "mmmmfmmmm", "mmmmfmfmm", "mmmmffmmm", "mmmmfffmm", 
                                             "mmmfmmmmm", "mmmfmmfmm", "mmmfmfmmm", "mmmfmffmm", "mmmffmmmm", "mmmffmfmm", "mmmfffmmm", "mmmffffmm", 
                                             "mmfmmmmmm", "mmfmmmfmm", "mmfmmfmmm", "mmfmmffmm", "mmfmfmmmm", "mmfmfmfmm", "mmfmffmmm", "mmfmfffmm", 
                                             "mmffmmmmm", "mmffmmfmm", "mmffmfmmm", "mmffmffmm", "mmfffmmmm", "mmfffmfmm", "mmffffmmm", "mmfffffmm", 
                                             "mfmmmmmmm", "mfmmmmfmm", "mfmmmfmmm", "mfmmmffmm", "mfmmfmmmm", "mfmmfmfmm", "mfmmffmmm", "mfmmfffmm", 
                                             "mfmfmmmmm", "mfmfmmfmm", "mfmfmfmmm", "mfmfmffmm", "mfmffmmmm", "mfmffmfmm", "mfmfffmmm", "mfmffffmm", 
                                             "mffmmmmmm", "mffmmmfmm", "mffmmfmmm", "mffmmffmm", "mffmfmmmm", "mffmfmfmm", "mffmffmmm", "mffmfffmm", 
                                             "mfffmmmmm", "mfffmmfmm", "mfffmfmmm", "mfffmffmm", "mffffmmmm", "mffffmfmm", "mfffffmmm", "mffffffmm", 
                                             "mmmmmmmfm", "mmmmmmffm", "mmmmmfmfm", "mmmmmfffm", "mmmmfmmfm", "mmmmfmffm", "mmmmffmfm", "mmmmffffm", 
                                             "mmmfmmmfm", "mmmfmmffm", "mmmfmfmfm", "mmmfmfffm", "mmmffmmfm", "mmmffmffm", "mmmfffmfm", "mmmfffffm", 
                                             "mmfmmmmfm", "mmfmmmffm", "mmfmmfmfm", "mmfmmfffm", "mmfmfmmfm", "mmfmfmffm", "mmfmffmfm", "mmfmffffm", 
                                             "mmffmmmfm", "mmffmmffm", "mmffmfmfm", "mmffmfffm", "mmfffmmfm", "mmfffmffm", "mmffffmfm", "mmffffffm", 
                                             "mfmmmmmfm", "mfmmmmffm", "mfmmmfmfm", "mfmmmfffm", "mfmmfmmfm", "mfmmfmffm", "mfmmffmfm", "mfmmffffm", 
                                             "mfmfmmmfm", "mfmfmmffm", "mfmfmfmfm", "mfmfmfffm", "mfmffmmfm", "mfmffmffm", "mfmfffmfm", "mfmfffffm", 
                                             "mffmmmmfm", "mffmmmffm", "mffmmfmfm", "mffmmfffm", "mffmfmmfm", "mffmfmffm", "mffmffmfm", "mffmffffm", 
                                             "mfffmmmfm", "mfffmmffm", "mfffmfmfm", "mfffmfffm", "mffffmmfm", "mffffmffm", "mfffffmfm", "mfffffffm", 
                                             "mmmmmmmmf", "mmmmmmfmf", "mmmmmfmmf", "mmmmmffmf", "mmmmfmmmf", "mmmmfmfmf", "mmmmffmmf", "mmmmfffmf", 
                                             "mmmfmmmmf", "mmmfmmfmf", "mmmfmfmmf", "mmmfmffmf", "mmmffmmmf", "mmmffmfmf", "mmmfffmmf", "mmmffffmf", 
                                             "mmfmmmmmf", "mmfmmmfmf", "mmfmmfmmf", "mmfmmffmf", "mmfmfmmmf", "mmfmfmfmf", "mmfmffmmf", "mmfmfffmf", 
                                             "mmffmmmmf", "mmffmmfmf", "mmffmfmmf", "mmffmffmf", "mmfffmmmf", "mmfffmfmf", "mmffffmmf", "mmfffffmf", 
                                             "mfmmmmmmf", "mfmmmmfmf", "mfmmmfmmf", "mfmmmffmf", "mfmmfmmmf", "mfmmfmfmf", "mfmmffmmf", "mfmmfffmf", 
                                             "mfmfmmmmf", "mfmfmmfmf", "mfmfmfmmf", "mfmfmffmf", "mfmffmmmf", "mfmffmfmf", "mfmfffmmf", "mfmffffmf", 
                                             "mffmmmmmf", "mffmmmfmf", "mffmmfmmf", "mffmmffmf", "mffmfmmmf", "mffmfmfmf", "mffmffmmf", "mffmfffmf", 
                                             "mfffmmmmf", "mfffmmfmf", "mfffmfmmf", "mfffmffmf", "mffffmmmf", "mffffmfmf", "mfffffmmf", "mffffffmf", 
                                             "mmmmmmmff", "mmmmmmfff", "mmmmmfmff", "mmmmmffff", "mmmmfmmff", "mmmmfmfff", "mmmmffmff", "mmmmfffff", 
                                             "mmmfmmmff", "mmmfmmfff", "mmmfmfmff", "mmmfmffff", "mmmffmmff", "mmmffmfff", "mmmfffmff", "mmmffffff", 
                                             "mmfmmmmff", "mmfmmmfff", "mmfmmfmff", "mmfmmffff", "mmfmfmmff", "mmfmfmfff", "mmfmffmff", "mmfmfffff", 
                                             "mmffmmmff", "mmffmmfff", "mmffmfmff", "mmffmffff", "mmfffmmff", "mmfffmfff", "mmffffmff", "mmfffffff", 
                                             "mfmmmmmff", "mfmmmmfff", "mfmmmfmff", "mfmmmffff", "mfmmfmmff", "mfmmfmfff", "mfmmffmff", "mfmmfffff", 
                                             "mfmfmmmff", "mfmfmmfff", "mfmfmfmff", "mfmfmffff", "mfmffmmff", "mfmffmfff", "mfmfffmff", "mfmffffff", 
                                             "mffmmmmff", "mffmmmfff", "mffmmfmff", "mffmmffff", "mffmfmmff", "mffmfmfff", "mffmffmff", "mffmfffff", 
                                             "mfffmmmff", "mfffmmfff", "mfffmfmff", "mfffmffff", "mffffmmff", "mffffmfff", "mfffffmff", "mffffffff", 
                                             "fmmmmmmmm", "fmmmmmfmm", "fmmmmfmmm", "fmmmmffmm", "fmmmfmmmm", "fmmmfmfmm", "fmmmffmmm", "fmmmfffmm", 
                                             "fmmfmmmmm", "fmmfmmfmm", "fmmfmfmmm", "fmmfmffmm", "fmmffmmmm", "fmmffmfmm", "fmmfffmmm", "fmmffffmm", 
                                             "fmfmmmmmm", "fmfmmmfmm", "fmfmmfmmm", "fmfmmffmm", "fmfmfmmmm", "fmfmfmfmm", "fmfmffmmm", "fmfmfffmm", 
                                             "fmffmmmmm", "fmffmmfmm", "fmffmfmmm", "fmffmffmm", "fmfffmmmm", "fmfffmfmm", "fmffffmmm", "fmfffffmm", 
                                             "ffmmmmmmm", "ffmmmmfmm", "ffmmmfmmm", "ffmmmffmm", "ffmmfmmmm", "ffmmfmfmm", "ffmmffmmm", "ffmmfffmm", 
                                             "ffmfmmmmm", "ffmfmmfmm", "ffmfmfmmm", "ffmfmffmm", "ffmffmmmm", "ffmffmfmm", "ffmfffmmm", "ffmffffmm", 
                                             "fffmmmmmm", "fffmmmfmm", "fffmmfmmm", "fffmmffmm", "fffmfmmmm", "fffmfmfmm", "fffmffmmm", "fffmfffmm", 
                                             "ffffmmmmm", "ffffmmfmm", "ffffmfmmm", "ffffmffmm", "fffffmmmm", "fffffmfmm", "ffffffmmm", "fffffffmm", 
                                             "fmmmmmmfm", "fmmmmmffm", "fmmmmfmfm", "fmmmmfffm", "fmmmfmmfm", "fmmmfmffm", "fmmmffmfm", "fmmmffffm", 
                                             "fmmfmmmfm", "fmmfmmffm", "fmmfmfmfm", "fmmfmfffm", "fmmffmmfm", "fmmffmffm", "fmmfffmfm", "fmmfffffm", 
                                             "fmfmmmmfm", "fmfmmmffm", "fmfmmfmfm", "fmfmmfffm", "fmfmfmmfm", "fmfmfmffm", "fmfmffmfm", "fmfmffffm", 
                                             "fmffmmmfm", "fmffmmffm", "fmffmfmfm", "fmffmfffm", "fmfffmmfm", "fmfffmffm", "fmffffmfm", "fmffffffm", 
                                             "ffmmmmmfm", "ffmmmmffm", "ffmmmfmfm", "ffmmmfffm", "ffmmfmmfm", "ffmmfmffm", "ffmmffmfm", "ffmmffffm", 
                                             "ffmfmmmfm", "ffmfmmffm", "ffmfmfmfm", "ffmfmfffm", "ffmffmmfm", "ffmffmffm", "ffmfffmfm", "ffmfffffm", 
                                             "fffmmmmfm", "fffmmmffm", "fffmmfmfm", "fffmmfffm", "fffmfmmfm", "fffmfmffm", "fffmffmfm", "fffmffffm", 
                                             "ffffmmmfm", "ffffmmffm", "ffffmfmfm", "ffffmfffm", "fffffmmfm", "fffffmffm", "ffffffmfm", "ffffffffm", 
                                             "fmmmmmmmf", "fmmmmmfmf", "fmmmmfmmf", "fmmmmffmf", "fmmmfmmmf", "fmmmfmfmf", "fmmmffmmf", "fmmmfffmf", 
                                             "fmmfmmmmf", "fmmfmmfmf", "fmmfmfmmf", "fmmfmffmf", "fmmffmmmf", "fmmffmfmf", "fmmfffmmf", "fmmffffmf", 
                                             "fmfmmmmmf", "fmfmmmfmf", "fmfmmfmmf", "fmfmmffmf", "fmfmfmmmf", "fmfmfmfmf", "fmfmffmmf", "fmfmfffmf", 
                                             "fmffmmmmf", "fmffmmfmf", "fmffmfmmf", "fmffmffmf", "fmfffmmmf", "fmfffmfmf", "fmffffmmf", "fmfffffmf", 
                                             "ffmmmmmmf", "ffmmmmfmf", "ffmmmfmmf", "ffmmmffmf", "ffmmfmmmf", "ffmmfmfmf", "ffmmffmmf", "ffmmfffmf", 
                                             "ffmfmmmmf", "ffmfmmfmf", "ffmfmfmmf", "ffmfmffmf", "ffmffmmmf", "ffmffmfmf", "ffmfffmmf", "ffmffffmf", 
                                             "fffmmmmmf", "fffmmmfmf", "fffmmfmmf", "fffmmffmf", "fffmfmmmf", "fffmfmfmf", "fffmffmmf", "fffmfffmf", 
                                             "ffffmmmmf", "ffffmmfmf", "ffffmfmmf", "ffffmffmf", "fffffmmmf", "fffffmfmf", "ffffffmmf", "fffffffmf", 
                                             "fmmmmmmff", "fmmmmmfff", "fmmmmfmff", "fmmmmffff", "fmmmfmmff", "fmmmfmfff", "fmmmffmff", "fmmmfffff", 
                                             "fmmfmmmff", "fmmfmmfff", "fmmfmfmff", "fmmfmffff", "fmmffmmff", "fmmffmfff", "fmmfffmff", "fmmffffff", 
                                             "fmfmmmmff", "fmfmmmfff", "fmfmmfmff", "fmfmmffff", "fmfmfmmff", "fmfmfmfff", "fmfmffmff", "fmfmfffff", 
                                             "fmffmmmff", "fmffmmfff", "fmffmfmff", "fmffmffff", "fmfffmmff", "fmfffmfff", "fmffffmff", "fmfffffff", 
                                             "ffmmmmmff", "ffmmmmfff", "ffmmmfmff", "ffmmmffff", "ffmmfmmff", "ffmmfmfff", "ffmmffmff", "ffmmfffff", 
                                             "ffmfmmmff", "ffmfmmfff", "ffmfmfmff", "ffmfmffff", "ffmffmmff", "ffmffmfff", "ffmfffmff", "ffmffffff", 
                                             "fffmmmmff", "fffmmmfff", "fffmmfmff", "fffmmffff", "fffmfmmff", "fffmfmfff", "fffmffmff", "fffmfffff", 
                                             "ffffmmmff", "ffffmmfff", "ffffmfmff", "ffffmffff", "fffffmmff", "fffffmfff", "ffffffmff", "fffffffff")]))
  })
  
  res$siblings <- ko(KidsOf = KidsOf, p = res$parents)
  res$siblings <- lapply(seq_along(res$siblings), function(i) res$siblings[[i]][res$siblings[[i]] %ni% pid[[i]]])

  g1 <- ko(KidsOf = KidsOf, p = res$gparents)
  res$unclesaunts <- lapply(seq_along(g1), function(i) g1[[i]][g1[[i]] %ni% res$parents[[i]]])
  
  #res$firstcousins <- ko(KidsOf = KidsOf, p = res$unclesaunts)
  
  g2 <- ko(KidsOf = KidsOf, p = res$ggparents)
  res$gunclesaunts <- lapply(seq_along(g2), function(i) g2[[i]][g2[[i]] %ni% res$gparents[[i]]])
  
  g3 <- ko(KidsOf = KidsOf, p = res$gggparents)
  res$ggunclesaunts <- lapply(seq_along(g3), function(i) g3[[i]][g3[[i]] %ni% res$ggparents[[i]]])
  
  g4 <- ko(KidsOf = KidsOf, p = res$ggggparents)
  res$gggunclesaunts <- lapply(seq_along(g4), function(i) g4[[i]][g4[[i]] %ni% res$gggparents[[i]]])
  
  g5 <- ko(KidsOf = KidsOf, p = res$gggggparents)
  res$ggggunclesaunts <- lapply(seq_along(g5), function(i) g5[[i]][g5[[i]] %ni% res$ggggparents[[i]]])
  
  g6 <- ko(KidsOf = KidsOf, p = res$ggggggparents)
  res$gggggunclesaunts <- lapply(seq_along(g6), function(i) g6[[i]][g6[[i]] %ni% res$gggggparents[[i]]])
  
  g7 <- ko(KidsOf = KidsOf, p = res$gggggggparents)
  res$ggggggunclesaunts <- lapply(seq_along(g7), function(i) g7[[i]][g7[[i]] %ni% res$ggggggparents[[i]]])
  
  g8 <- ko(KidsOf = KidsOf, p = res$ggggggggparents)
  res$gggggggunclesaunts <- lapply(seq_along(g8), function(i) g8[[i]][g8[[i]] %ni% res$gggggggparents[[i]]])
  
  # res$spouse <- so(opop = opop, p = pid)
  # 
  # res$children <- ko(KidsOf = KidsOf, p = pid)
 
  for (i in 1:length(names(res))) {
    res[[i]][sapply(res[[i]], function(x) length(x) == 0)] <- NA
  }

  return(res)
}