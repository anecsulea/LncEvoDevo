########################################################################

files=c("resample.reads.pc.lncRNAs.R", "differential.expression.ages.R", "differential.expression.ages.resampled.pc.lncRNA.R", "differential.expression.consecutive.stages.pc.lnc.R")  

########################################################################

for(file in files){
  print(paste("running ",file))
  
  load=T
  process=T
  load.annot=T
  load.coverage=T
  
  source(file)
  rm(list=setdiff(ls(), c("files", "file")))
}

########################################################################
