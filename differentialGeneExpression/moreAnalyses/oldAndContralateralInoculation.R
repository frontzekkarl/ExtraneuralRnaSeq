# comparison of main cohort to mice inoculated at 1 year of age, "PBO" and to young mice
# inoculated in the contralateral hemisphere "PBH"

# read and prepare blood 8 wpi and muscle 8 wpi from main cohort according to:
# analysisBloodMainCohort.R and analysisMuscleMainCohort.R from differentialGeneExpression/ 

library(readxl)

####### PBO

PBO_mu <- read.delim("Muscle/result--RML6_old_8--over--NBH_old_8.txt") # PBO edgeR output
PBO_mu <- subset(PBO_mu,(PBO_mu$fdr<0.05) & (abs(PBO_mu$log2.Ratio)>0.5),c(3,13,14,21))

# intersect PBO DEGs with main cohort DEGs
PBO_mu$gene_name %in% mu_8$gene_name # 0


PB0_bl <- read_xlsx("Blood/result--RML6_old_8--over--NBH_old_8.xlsx")
PB0_bl <- subset(PB0_bl,(PB0_bl$fdr<0.05) & (abs(PB0_bl$`log2 Ratio`>0.5)))
# 0 DEG in Blood from PBO
PB0_bl[(PB0_bl$'fdr' < 0.05) & (abs(PB0_bl$`log2 Ratio`>0.5)),] # 0       

### PBH_mu <- read.delim("Muscle/CohortPBH/result--RML6_8--over--NBH_8.txt")
PBH_mu <- read.delim("Muscle/CohortPBH/result--RML6_8--over--NBH_8.txt")
PBH_mu <- subset(PBH_mu,(PBH_mu$fdr<0.05) & (abs(PBH_mu$log2.Ratio>0.5)),c(3,13,14,21))
PBH_mu[(PBH_mu$fdr < 0.05) & (abs(PBH_mu$log2.Ratio>0.5)),] # 7
PBH_mu$gene_name %in% mu_8$gene_name  #0

PBH_bl <- read_xlsx("Blood/CohortPBH/result--RML6_8--over--NBH_8.xlsx")
PBH_bl <- subset(PBH_bl,(fdr<0.05)& (abs(PBH_bl$`log2 Ratio`)>0.5))
PBH_bl # 0
