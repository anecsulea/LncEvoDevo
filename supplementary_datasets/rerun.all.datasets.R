########################################################################

all.files=system("ls | grep ^data | grep R$", intern=T)

########################################################################

for(file in all.files){
  print(file)
  
 
  source(file)
  rm(list=setdiff(ls(), c("files", "file")))
  
}

########################################################################
