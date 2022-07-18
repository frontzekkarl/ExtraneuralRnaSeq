setwd("./rawData/mainCohort/blood/")
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(patchwork)

bl_4 <- read.delim("result--RML6_4--over--NBH_4.txt")
bl_8 <- read.delim("result--RML6_8--over--NBH_8.txt")
bl_12 <- read.delim("result--RML6_12--over--NBH_12.txt")
bl_14 <- read.delim("result--RML6_14--over--NBH_14.txt")
bl_16 <- read.delim("result--RML6_16--over--NBH_16.txt")
bl_18 <- read.delim("result--RML6_18--over--NBH_18.txt")
bl_20 <- read.delim("result--RML6_20--over--NBH_20.txt")
bl_term <- read.delim("result--RML6_term--over--NBH_term.txt")
#bl_none <- read.delim("result--none_36wks_---over--none_7wks_-.txt")

# test number of DE genes 
bl_4_summary <- summary(bl_4$fdr<0.05 & abs(bl_4$log2.Ratio)>0.5) 
print(paste0('number of DE genes in blood 4 weeks: ',bl_4_summary['TRUE'])) # 319
bl_8_summary <- summary(bl_8$fdr<0.05 & abs(bl_8$log2.Ratio)>0.5) 
print(paste0('number of DE genes in blood 8 weeks: ',bl_8_summary['TRUE'])) # 0
bl_12_summary <- summary(bl_12$fdr<0.05 & abs(bl_12$log2.Ratio)>0.5) 
print(paste0('number of DE genes in blood 12 weeks: ',bl_12_summary['TRUE'])) # 0
bl_14_summary <- summary(bl_14$fdr<0.05 & abs(bl_14$log2.Ratio)>0.5) 
print(paste0('number of DE genes in blood 14 weeks: ',bl_14_summary['TRUE'])) # 100
bl_16_summary <- summary(bl_16$fdr<0.05 & abs(bl_16$log2.Ratio)>0.5) 
print(paste0('number of DE genes in blood 16 weeks: ',bl_16_summary['TRUE'])) # 5
bl_18_summary <- summary(bl_18$fdr<0.05 & abs(bl_18$log2.Ratio)>0.5) 
print(paste0('number of DE genes in blood 18 weeks: ',bl_18_summary['TRUE'])) # 201
bl_20_summary <- summary(bl_20$fdr<0.05 & abs(bl_20$log2.Ratio)>0.5) 
print(paste0('number of DE genes in blood 20 weeks: ',bl_20_summary['TRUE'])) # 8
bl_term_summary <- summary(bl_term$fdr<0.05 & abs(bl_term$log2.Ratio)>0.5) 
print(paste0('number of DE genes in blood term weeks: ',bl_term_summary['TRUE'])) # 16

# subset DEG
bl_4_genes <- subset(bl_4,(bl_4$fdr<0.05) & (abs(bl_4$log2.Ratio)>0.5),c(3,13,14,21))
bl_14_genes <- subset(bl_14,(bl_14$fdr<0.05) & (abs(bl_14$log2.Ratio)>0.5),c(3,13,14,21))
bl_16_genes <- subset(bl_16,(bl_16$fdr<0.05) & (abs(bl_16$log2.Ratio)>0.5),c(3,13,14,21))
bl_18_genes <- subset(bl_18,(bl_18$fdr<0.05) & (abs(bl_18$log2.Ratio)>0.5),c(3,13,14,21))
bl_20_genes <- subset(bl_20,(bl_20$fdr<0.05) & (abs(bl_20$log2.Ratio)>0.5),c(3,13,14,21))
bl_term_genes <- subset(bl_term,(bl_term$fdr<0.05) & (abs(bl_term$log2.Ratio)>0.5),c(3,13,14,21))

paste0('number of upregulated genes week 4: ',length(bl_4_genes$log2.Ratio[bl_4_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 4: ',length(bl_4_genes$log2.Ratio[bl_4_genes$log2.Ratio<0]))

paste0('number of upregulated genes week 14: ',length(bl_14_genes$log2.Ratio[bl_14_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 14: ',length(bl_14_genes$log2.Ratio[bl_14_genes$log2.Ratio<0]))

paste0('number of upregulated genes week 16: ',length(bl_16_genes$log2.Ratio[bl_16_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 16: ',length(bl_16_genes$log2.Ratio[bl_16_genes$log2.Ratio<0]))

paste0('number of upregulated genes week 18: ',length(bl_18_genes$log2.Ratio[bl_18_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 18: ',length(bl_18_genes$log2.Ratio[bl_18_genes$log2.Ratio<0]))

paste0('number of upregulated genes week 20: ',length(bl_20_genes$log2.Ratio[bl_20_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 20: ',length(bl_20_genes$log2.Ratio[bl_20_genes$log2.Ratio<0]))

paste0('number of upregulated genes week term: ',length(bl_term_genes$log2.Ratio[bl_term_genes$log2.Ratio>0]))
paste0('number of downregulated genes week term: ',length(bl_term_genes$log2.Ratio[bl_term_genes$log2.Ratio<0]))





# show genes
# columns 1,2,3,12,13,14,21

# week 4
rowind_4 <- which(!grepl("NA", rownames(bl_4[(bl_4$fdr<0.05 & abs(bl_4$log2.Ratio)>0.5),])),arr.ind=TRUE) 
bl_4_genes <- bl_4[rowind_4,c(3,13,14,21)]
# week 14
rowind_14 <-  which(!grepl("NA", rownames(bl_14[(bl_14$fdr<0.05 & abs(bl_14$log2.Ratio)>0.5),])),arr.ind=TRUE) 
bl_14_genes <- bl_14[rowind_14,c(3,13,14,21)]
# week 16
rowind_16 <-  which(!grepl("NA", rownames(bl_16[(bl_16$fdr<0.05 & abs(bl_16$log2.Ratio)>0.5),])),arr.ind=TRUE) 
bl_16_genes <- bl_16[rowind_16,c(3,13,14,21)]
# week 18
rowind_18 <-  which(!grepl("NA", rownames(bl_18[(bl_18$fdr<0.05 & abs(bl_18$log2.Ratio)>0.5),])),arr.ind=TRUE) 
bl_18_genes <- bl_18[rowind_18,c(3,13,14,21)]
# week 20
rowind_20 <-  which(!grepl("NA", rownames(bl_20[(bl_20$fdr<0.05 & abs(bl_20$log2.Ratio)>0.5),])),arr.ind=TRUE) 
bl_20_genes <- bl_20[rowind_20,c(3,13,14,21)]
# term
rowind_term <-  which(!grepl("NA", rownames(bl_term[(bl_term$fdr<0.05 & abs(bl_term$log2.Ratio)>0.5),])),arr.ind=TRUE) 
bl_term_genes <- bl_term[rowind_term,c(3,13,14,21)]

# DEG up&down

sum(bl_4_genes$log2.Ratio>0)
sum(bl_4_genes$log2.Ratio<0)

sum(bl_14_genes$log2.Ratio>0)
sum(bl_14_genes$log2.Ratio<0)

sum(bl_16_genes$log2.Ratio>0)
sum(bl_16_genes$log2.Ratio<0)

sum(bl_18_genes$log2.Ratio>0)
sum(bl_18_genes$log2.Ratio<0)

sum(bl_20_genes$log2.Ratio>0)
sum(bl_20_genes$log2.Ratio<0)

sum(bl_term_genes$log2.Ratio>0)
sum(bl_term_genes$log2.Ratio<0)

# write csv
write.csv2(bl_4_genes,"Analysis/bl4_genes_blood.csv")
write.csv2(bl_14_genes,"Analysis/bl14_genes_blood.csv")
write.csv2(bl_16_genes,"Analysis/bl16_genes_blood.csv")
write.csv2(bl_18_genes,"Analysis/bl18_genes_blood.csv")
write.csv2(bl_20_genes,"Analysis/bl20_genes_blood.csv")
write.csv2(bl_term_genes,"Analysis/blterm_genes_blood.csv")
