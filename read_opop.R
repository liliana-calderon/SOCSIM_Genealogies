## Function written by Diego Alburez-Gutierrez. 
# This will be later integrated in the rsocsim package

read_opop <- function(path){
  opop <- read.table(file=path,header=F,as.is=T)  
  ## assign names to columns
  names(opop)<-c("pid","fem","group",
                 "nev","dob","mom","pop","nesibm","nesibp",
                 "lborn","marid","mstat","dod","fmult")
  return(opop)
}