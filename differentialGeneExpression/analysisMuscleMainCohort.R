# identify differentially expressed genes in muscle from main cohort

library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(patchwork)
mu_4 <- read.delim("result--RML6_4--over--NBH_4.txt")
mu_8 <- read.delim("result--RML6_8--over--NBH_8.txt")
mu_12 <- read.delim("result--RML6_12--over--NBH_12.txt")
mu_14 <- read.delim("result--RML6_14--over--NBH_14.txt")
mu_16 <- read.delim("result--RML6_16--over--NBH_16.txt")
mu_18 <- read.delim("result--RML6_18--over--NBH_18.txt")
mu_20 <- read.delim("result--RML6_20--over--NBH_20.txt")
mu_term <- read.delim("result--RML6_term--over--NBH_term.txt")



mu_4_summary <- summary(mu_4$fdr<0.05 & abs(mu_4$log2.Ratio)>0.5) 
print(paste0('number of DE genes in muood 4 weeks: ',mu_4_summary['TRUE'])) # 26
mu_8_summary <- summary(mu_8$fdr<0.05 & abs(mu_8$log2.Ratio)>0.5) 
print(paste0('number of DE genes in muood 8 weeks: ',mu_8_summary['TRUE'])) # 106
mu_12_summary <- summary((mu_12$fdr<0.05) & (abs(mu_12$log2.Ratio)>0.5)) 
print(paste0('number of DE genes in muood 12 weeks: ',mu_12_summary['TRUE'])) # 55
mu_14_summary <- summary(mu_14$fdr<0.05 & abs(mu_14$log2.Ratio)>0.5) 
print(paste0('number of DE genes in muood 14 weeks: ',mu_14_summary['TRUE'])) # 476
mu_16_summary <- summary(mu_16$fdr<0.05 & abs(mu_16$log2.Ratio)>0.5) 
print(paste0('number of DE genes in muood 16 weeks: ',mu_16_summary['TRUE'])) # 0
mu_18_summary <- summary(mu_18$fdr<0.05 & abs(mu_18$log2.Ratio)>0.5) 
print(paste0('number of DE genes in muood 18 weeks: ',mu_18_summary['TRUE'])) # 291
mu_20_summary <- summary(mu_20$fdr<0.05 & abs(mu_20$log2.Ratio)>0.5) 
print(paste0('number of DE genes in muood 20 weeks: ',mu_20_summary['TRUE'])) # 197
mu_term_summary <- summary(mu_term$fdr<0.05 & abs(mu_term$log2.Ratio)>0.5) 
print(paste0('number of DE genes in muood term weeks: ',mu_term_summary['TRUE'])) # 1320

mu_4_genes <- subset(mu_4,(mu_4$fdr<0.05) & (abs(mu_4$log2.Ratio)>0.5),c(3,13,14,21))
mu_8_genes <- subset(mu_8,(mu_8$fdr<0.05) & (abs(mu_8$log2.Ratio)>0.5),c(3,13,14,21))
mu_12_genes <- subset(mu_12,(mu_12$fdr<0.05) & (abs(mu_12$log2.Ratio)>0.5),c(3,13,14,21))
mu_14_genes <- subset(mu_14,(mu_14$fdr<0.05) & (abs(mu_14$log2.Ratio)>0.5),c(3,13,14,21))
mu_18_genes <- subset(mu_18,(mu_18$fdr<0.05) & (abs(mu_18$log2.Ratio)>0.5),c(3,13,14,21))
mu_20_genes <- subset(mu_20,(mu_20$fdr<0.05) & (abs(mu_20$log2.Ratio)>0.5),c(3,13,14,21))
mu_term_genes <- subset(mu_term,(mu_term$fdr<0.05) & (abs(mu_term$log2.Ratio)>0.5),c(3,13,14,21))
paste0('number of upregulated genes week 4: ',length(mu_4_genes$log2.Ratio[mu_4_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 4: ',length(mu_4_genes$log2.Ratio[mu_4_genes$log2.Ratio<0]))

paste0('number of upregulated genes week 8: ',length(mu_8_genes$log2.Ratio[mu_8_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 8: ',length(mu_8_genes$log2.Ratio[mu_8_genes$log2.Ratio<0]))

paste0('number of upregulated genes week 12: ',length(mu_12_genes$log2.Ratio[mu_12_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 12: ',length(mu_12_genes$log2.Ratio[mu_12_genes$log2.Ratio<0]))

paste0('number of upregulated genes week 14: ',length(mu_14_genes$log2.Ratio[mu_14_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 14: ',length(mu_14_genes$log2.Ratio[mu_14_genes$log2.Ratio<0]))

paste0('number of upregulated genes week 18: ',length(mu_18_genes$log2.Ratio[mu_18_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 18: ',length(mu_18_genes$log2.Ratio[mu_18_genes$log2.Ratio<0]))

paste0('number of upregulated genes week 20: ',length(mu_20_genes$log2.Ratio[mu_20_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 20: ',length(mu_20_genes$log2.Ratio[mu_20_genes$log2.Ratio<0]))

paste0('number of upregulated genes week term: ',length(mu_term_genes$log2.Ratio[mu_term_genes$log2.Ratio>0]))
paste0('number of downregulated genes week term: ',length(mu_term_genes$log2.Ratio[mu_term_genes$log2.Ratio<0]))


# show genes
# columns 1,2,3,12,13,14,21

# week 4
rowind_4 <- which(!grepl("NA", rownames(mu_4[(mu_4$fdr<0.05 & abs(mu_4$log2.Ratio)>0.5),])),arr.ind=TRUE) 
mu_4_genes <- mu_4[rowind_4,c(3,13,14,21)]
# week 8
rowind_8 <- which(!grepl("NA", rownames(mu_8[(mu_8$fdr<0.05 & abs(mu_8$log2.Ratio)>0.5),])),arr.ind=TRUE) 
mu_8_genes <- mu_8[rowind_8,c(3,13,14,21)]
# week 12
rowind_12 <- which(!grepl("NA", rownames(mu_12[(mu_12$fdr<0.05 & abs(mu_12$log2.Ratio)>0.5),])),arr.ind=TRUE) 
mu_12_genes <- mu_12[rowind_12,c(3,13,14,21)]
# week 14
rowind_14 <-  which(!grepl("NA", rownames(mu_14[(mu_14$fdr<0.05 & abs(mu_14$log2.Ratio)>0.5),])),arr.ind=TRUE) 
mu_14_genes <- mu_14[rowind_14,c(3,13,14,21)]
# week 18
rowind_18 <-  which(!grepl("NA", rownames(mu_18[(mu_18$fdr<0.05 & abs(mu_18$log2.Ratio)>0.5),])),arr.ind=TRUE) 
mu_18_genes <- mu_18[rowind_18,c(3,13,14,21)]
# week 20
rowind_20 <-  which(!grepl("NA", rownames(mu_20[(mu_20$fdr<0.05 & abs(mu_20$log2.Ratio)>0.5),])),arr.ind=TRUE) 
mu_20_genes <- mu_20[rowind_20,c(3,13,14,21)]
# term
rowind_term <-  which(!grepl("NA", rownames(mu_term[(mu_term$fdr<0.05 & abs(mu_term$log2.Ratio)>0.5),])),arr.ind=TRUE) 
mu_term_genes <- mu_term[rowind_term,c(3,13,14,21)]
mu_20_genes

# write csv
write.csv2(mu_4_genes,"Analysis/mu4_genes_muscle.csv")
write.csv2(mu_8_genes,"Analysis/mu8_genes_muscle.csv")
write.csv2(mu_12_genes,"Analysis/mu12_genes_muscle.csv")
write.csv2(mu_14_genes,"Analysis/mu14_genes_muscle.csv")
write.csv2(mu_18_genes,"Analysis/mu18_genes_muscle.csv")
write.csv2(mu_20_genes,"Analysis/mu20_genes_muscle.csv")
write.csv2(mu_term_genes,"Analysis/muterm_genes_muscle.csv")
summary(abs(mu_12_genes$log2.Ratio) < 0.5)

sum(mu_4_genes$log2.Ratio>0)
sum(mu_4_genes$log2.Ratio<0)

sum(mu_8_genes$log2.Ratio>0)
sum(mu_8_genes$log2.Ratio<0)

sum(mu_12_genes$log2.Ratio>0)
sum(mu_12_genes$log2.Ratio<0)

sum(mu_14_genes$log2.Ratio>0)
sum(mu_14_genes$log2.Ratio<0)

sum(mu_18_genes$log2.Ratio>0)
sum(mu_18_genes$log2.Ratio<0)

sum(mu_20_genes$log2.Ratio>0)
sum(mu_20_genes$log2.Ratio<0)

sum(mu_term_genes$log2.Ratio>0)
sum(mu_term_genes$log2.Ratio<0)
