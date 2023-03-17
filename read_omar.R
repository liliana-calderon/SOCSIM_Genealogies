## Function written by Diego Alburez-Gutierrez. 
# This will be later integrated in the rsocsim package

read_omar <- function(path){
  omar<-read.table(file = path, header = F, as.is = T)
  names(omar)<-c("mid","wpid","hpid","dstart","dend", "rend","wprior","hprior")
  return(omar)
}