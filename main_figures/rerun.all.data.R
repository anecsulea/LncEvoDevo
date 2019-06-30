########################################################################

all.files=system("ls | grep ^data | grep R$", intern=T)

########################################################################

for(file in all.files){
  print(file)
  
  if(file=="data.expression.conservation.R"){
    for(method in c("pearson", "spearman")){
      source(file)
      rm(list=setdiff(ls(), c("files", "file", "method")))
    }
  } else{
    source(file)
    rm(list=setdiff(ls(), c("files", "file")))
  }
}

########################################################################
