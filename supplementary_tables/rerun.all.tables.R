########################################################################

all.files=c("table.cell.markers.R", "table.nbexplnc.R", "table.all.conservation.R", "table.ortho.mrc.R", "table.sequence.conservation.R")

########################################################################

for(file in all.files){
  print(file)
  
 
  source(file)
  rm(list=setdiff(ls(), c("files", "file")))
  
}

########################################################################
