# differential editing analysis from high confidence RNA editing sites 

source("/home/ubuntu/Bulk_RNAseq/RediEli/REDIT_LLR.R") # read Redit
library(tidyverse)
bloodRedi <- read.delim("/home/ubuntu/Bulk_RNAseq/RediEli/brainRediEli.csv",sep=",")
#bloodRedi <- as_tibble(bloodRedi)
bloodRedi
wpiPosUnique <- unique(bloodRedi$wpiPosition) # select unique editing sites
dF <- t(subset(bloodRedi,wpiPosition==wpiPosUnique[1], select=c(G,A)))
dF
#subset(bloodRedi,wpiPosition=="term+172092348", select=Treatment) 

dP<- data.frame(site=c(),pValue=c()) #preassign empty dataframe

for (i in 1:length(wpiPosUnique)) {
  wpiPos <- wpiPosUnique[i] 
  dF <- t(subset(bloodRedi,wpiPosition==wpiPos, select=c(G,A))) # generate table of G and A numbers
  treat <- subset(bloodRedi,wpiPosition==wpiPos, select=Treatment_y) # generate treatment vector
  #if(wpiPos=="term+41654252") next
  #print(wpiPos)
  redit <- REDIT_LLR(data=dF, groups = as.vector(unlist(treat))) # redit test
  #print(length(treat$Treatment))
  dT <- data.frame(site=wpiPos,pValue=redit$p.value)
  dP <- rbind(dP,dT) # append data
}

write_csv(dP,'/home/ubuntu/Bulk_RNAseq/RediEli/brainEditingSitesRedit.csv')
