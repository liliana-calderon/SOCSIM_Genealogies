#----------------------------------------------------------------------------------------------------
# Function to retrieve members of a kin network of a given ego from SOCSIM microsimulation outputs
# considering direct ancestors up to 8th degree of consanguinity and their offspring

# To run the functions, .opop and omar files must be set in the GlobalEnv

# This is a modified version of the retrieve_kin function created by Mallika Snyder, included in the rsocsim package. 
# Here we add the direct ancestors from 4th to 8th degree of consanguinity and their offspring
# The name of the vectors for grandparents and great-grand parents is also different, 

# Created by Liliana Calderon on 23-09-2022
# Last modified by Liliana Calderon on 31-07-2023
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


  
  res$siblings <- ko(KidsOf = KidsOf, p = res$parents)
  res$siblings <- lapply(seq_along(res$siblings), function(i) res$siblings[[i]][res$siblings[[i]] %ni% pid[[i]]])

  g1 <- ko(KidsOf = KidsOf, p = res$gparents)
  res$unclesaunts <- lapply(seq_along(g1), function(i) g1[[i]][g1[[i]] %ni% res$parents[[i]]])
  
  res$firstcousins <- ko(KidsOf = KidsOf, p = res$unclesaunts)
  
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
  
  # res$spouse <- so(opop = opop, p = pid)
  # 
  # res$children <- ko(KidsOf = KidsOf, p = pid)
 
  for (i in 1:length(names(res))) {
    res[[i]][sapply(res[[i]], function(x) length(x) == 0)] <- NA
  }

  return(res)
}