########################################################################

all.files=system("ls | grep figure | grep R$", intern=T)
all.files=setdiff(all.files, "rerun.all.figures.R")

########################################################################

for(file in all.files){
  print(file)
  
  source(file)
  rm(list=setdiff(ls(), c("files", "file")))
  
}

########################################################################
