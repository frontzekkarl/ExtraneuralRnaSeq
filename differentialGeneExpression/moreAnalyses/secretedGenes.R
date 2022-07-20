### intersect blood DEGs with secreted genes from:
## The Metazoa (Human/Animal) Secretome and Subcellular Proteome KnowledgeBase (MetazSecKB, (Meinken et al., 2015))
# here imported as "funseckb2_search_results.txt"

sec <- read.csv("funseckb2_search_results.txt",header=FALSE,sep="\t")
sec <- as.character(sec$V1)

## convert uniprot -> gene id
library(org.Mm.eg.db)
columns(org.Mm.eg.db)
sec_symbol <- select(org.Mm.eg.db, sec, "SYMBOL", "UNIPROT")
sec <- sec_symbol$SYMBOL
sec <- sec[!is.na(sec)] # remove na

## load blood genes
###################

setwd("/home/ubuntu/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/BloodPlusBatch")
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

bl_term_genes
## test intersect
sec_4 <- bl_4_genes[bl_4_genes$gene_name %in% sec,]
sec_4 <- sec_4[order(-sec_4$log2.Ratio),]
sec_14 <- bl_14_genes[bl_14_genes$gene_name %in% sec,]
sec_14 <- sec_14[order(-sec_14$log2.Ratio),]
sec_16 <- bl_16_genes[bl_16_genes$gene_name %in% sec,]
sec_16 <- sec_16[order(-sec_16$log2.Ratio),]
sec_18 <- bl_18_genes[bl_18_genes$gene_name %in% sec,]
sec_18 <- sec_18[order(-sec_18$log2.Ratio),]
sec_20 <- bl_20_genes[bl_20_genes$gene_name %in% sec,]
sec_20 <- sec_20[order(-sec_20$log2.Ratio),]
sec_term <- bl_term_genes[bl_term_genes$gene_name %in% sec,]
sec_term <- sec_term[order(-sec_term$log2.Ratio),]

sec_4$timepoint <- 4
sec_14$timepoint <- 14
sec_16$timepoint <- 16
sec_18$timepoint <- 18
sec_20$timepoint <- 20
sec_term$timepoint <- 22

# write secreted genes
sec_bl <- do.call(rbind, list(sec_4,sec_14,sec_16,sec_18,sec_20,sec_term))
setwd("/home/ubuntu/Bulk_RNAseq/")
write.table(sec_bl,"secreted_genes_blood.csv",sep="\t",quote=FALSE)

sec_4$gene_name <- as.character(sec_4$gene_name)
sec_14$gene_name <- as.character(sec_14$gene_name)
sec_16$gene_name <- as.character(sec_16$gene_name)
sec_18$gene_name <- as.character(sec_18$gene_name)
sec_20$gene_name <- as.character(sec_20$gene_name)
sec_term$gene_name <- as.character(sec_term$gene_name)


# count monotonic genes
sec_blood <- table(c(sec_4$gene_name,sec_14$gene_name,sec_16$gene_name,sec_18$gene_name,sec_20$gene_name,sec_term$gene_name))
sec_blood[sec_blood > 1]

Reduce(intersect,list(c(sec_4$gene_name,sec_14$gene_name,sec_16$gene_name,sec_18$gene_name,sec_20$gene_name,sec_term$gene_name)))

c(sec_4$gene_name,sec_14$gene_name,sec_16$gene_name,sec_18$gene_name,sec_20$gene_name,sec_term$gene_name)

