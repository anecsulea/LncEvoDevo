########################################################

path="LncEvoDevo/"
pathFigures=paste(path, "figures/", sep="")
pathDatasets=paste(path, "supplementary_datasets/", sep="")
pathTables=paste(path, "supplementary_tables/", sep="")

########################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#########################################################

col.tissues=c("dodgerblue4", "seagreen", "indianred4","goldenrod")
names(col.tissues)=c("Brain", "Kidney", "Liver","Testis")

#########################################################

tissue.order=c("Brain", "Kidney", "Liver", "Testis")
age.order=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")
age4.order=c("EarlyEmbryo", "LateEmbryo", "Newborn", "AdultAged")

shortname.age=c("mid-stage embryo", "late embryo", "newborn", "young adult", "aged adult")
names(shortname.age)=age.order

#########################################################

tissage.order=kronecker(tissue.order, age.order, paste, sep="_")
shortname.tiss=c("brain", "kidney", "liver", "testes")

names(shortname.tiss)=tissue.order

#########################################################


## col.species=c("red", "navy", "seagreen")
## names(col.species)=c("Mouse", "Rat", "Chicken")

#########################################################

col.ages=c(rgb(216,216,216, maxColorValue=255), rgb(165, 165, 165, maxColorValue=255), rgb(126, 126, 126, maxColorValue=255), rgb(66, 66, 66, maxColorValue=255), rgb(16, 16, 16, maxColorValue=255)) 
names(col.ages)=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")


syn.age=c("mid-stage em.", "late embryo", "newborn", "young adult", "aged adult")
longsyn.age=c("mid-stage embryo", "late embryo", "newborn", "young adult", "aged adult")

#########################################################

pch.smallage=c(21:24)
names(pch.smallage)=c("EarlyEmbryo", "LateEmbryo", "Newborn", "AdultAged")

pch.age=c(21:25)
names(pch.age)=c("EarlyEmbryo", "LateEmbryo", "Newborn", "Adult", "Aged")

#########################################################

col.tissage=c("cadetblue1", "steelblue1", "dodgerblue3", "blue1", "blue4",  "aquamarine", "seagreen2",  "springgreen4", "chartreuse4", "darkgreen", "lightsalmon",  "tomato1", "red",  "brown",  "red4", "yellow", "gold", "goldenrod", "darkorange4")

names(col.tissage)=setdiff(kronecker(tissue.order, age.order, paste, sep="_"), "Testis_EarlyEmbryo")

#########################################################

pch.allsp=c(21, 23, 24)
names(pch.allsp)=c("mouse", "rat", "chicken")

col.allsp=c("black", "gray40", "gray80")
names(col.allsp)=c("mouse", "rat", "chicken")

#########################################################
