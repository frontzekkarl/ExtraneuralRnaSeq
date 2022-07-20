library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(patchwork)
sp_4 <- read.delim("Spleen/result--RML6_4--over--NBH_4.txt")
sp_8 <- read.delim("Spleen/result--RML6_8--over--NBH_8.txt")
sp_12 <- read.delim("Spleen/result--RML6_12--over--NBH_12.txt")
sp_14 <- read.delim("Spleen/result--RML6_14--over--NBH_14.txt")
sp_16 <- read.delim("Spleen/result--RML6_16--over--NBH_16.txt")
sp_18 <- read.delim("Spleen/result--RML6_18--over--NBH_18.txt")
sp_20 <- read.delim("Spleen/result--RML6_20--over--NBH_20.txt")
sp_term <- read.delim("Spleen/result--RML6_term--over--NBH_term.txt")
sp_none <- read.delim("Spleen/result--none_36wks_---over--none_7wks_-.txt")


# number up/down
#  weeks 8,12,14,18,20 don't show deg, here already excluded

paste0('number of upregulated genes week 4: ',length(sp_4_genes$log2.Ratio[sp_4_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 4: ',length(sp_4_genes$log2.Ratio[sp_4_genes$log2.Ratio<0.]))

paste0('number of upregulated genes week 16: ',length(sp_16_genes$log2.Ratio[sp_16_genes$log2.Ratio>0]))
paste0('number of downregulated genes week 16: ',length(sp_16_genes$log2.Ratio[sp_16_genes$log2.Ratio<0]))

paste0('number of upregulated genes week term: ',length(sp_term_genes$log2.Ratio[sp_term_genes$log2.Ratio>0]))
paste0('number of downregulated genes week term: ',length(sp_term_genes$log2.Ratio[sp_term_genes$log2.Ratio<0]))


#write genes
write.csv2(sp_4_genes,"Spleen/Analysis/sp4_genes.csv")
write.csv2(sp_16_genes,"Spleen/Analysis/sp_16_genes.csv")
write.csv2(sp_term_genes,"Spleen/Analysis/term_genes.csv")



# show genes
# columns 1,2,3,12,13,14,21
# week 4

rowind_4 <- which(!grepl("NA", rownames(sp_4[(sp_4$fdr<0.05 & abs(sp_4$log2.Ratio)>0.5),])),arr.ind=TRUE) 
sp_4_genes <- sp_4[rowind_4,c(3,13,14,21)]

# week 16
rowind_16 <-  which(!grepl("NA", rownames(sp_16[(sp_16$fdr<0.05 & abs(sp_16$log2.Ratio)>0.5),])),arr.ind=TRUE) 
sp_16_genes <- sp_16[rowind_16,c(3,13,14,21)]

# order by log2 fold change
sp_16_genes[order(sp_16_genes$log2.Ratio),]

# terminal
rowind_term <-  which(!grepl("NA", rownames(sp_term[(sp_term$fdr<0.05 & abs(sp_term$log2.Ratio)>0.5),])),arr.ind=TRUE) 
sp_term_genes <- sp_term[rowind_term,c(3,13,14,21)]

# # of up- and downregulated DEG
sum(sp_16_genes$log2.Ratio>0)
sum(sp_16_genes$log2.Ratio<0)

sum(sp_term_genes$log2.Ratio>0)
sum(sp_term_genes$log2.Ratio<0)
# write csv

# make list of genes and log2fc
# term
#--------------

term2 <- subset(sp_term_genes[order(sp_term_genes$log2.Ratio),], select = c(gene_name,log2.Ratio))
sp16_2 <- subset(sp_16_genes[order(sp_16_genes$log2.Ratio),], select = c(gene_name,log2.Ratio))
sp4_2 <- subset(sp_4_genes[order(sp_16_genes$log2.Ratio),], select = c(gene_name,log2.Ratio))
