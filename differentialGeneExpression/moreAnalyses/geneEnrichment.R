# perform gene enrichment analysis and gene ontology annotation analysis on blood 4 wpi samples

library(clusterProfiler)
library(dplyr)
library(enrichplot)

### bl4
## -----
#load data
bl_4 <- read.delim("result--RML6_4--over--NBH_4.txt") # corresponds to edgeR output from blood 4 wpi

# get significantly enriched genes
bl_4_summary <- summary(bl_4$fdr<0.05 & abs(bl_4$log2.Ratio)>0.5) 

print(paste0('number of DE genes in blood 4 weeks: ',bl_4_summary['TRUE'])) # 319

# make list of all genes and enriched genes
bl_4_gene_list <- subset(bl_4,(bl_4$fdr<0.05) & (abs(bl_4$log2.Ratio)>0.5),c(3,13,14,21))
bl4_short <- bl_4[,c(3,13,14,20,21)]
bl4_sig_genelist <- bl_4_gene_list$gene_name
bl4_all_genelist <- subset(bl4_short,bl4_short$isPresent == TRUE, select=gene_name)
length(bl4_all_genelist$gene_name)

# perform enrichGO on gene names
bl4_enrichGO <- enrichGO(gene          = bl4_sig_genelist,
                          universe      = bl4_all_genelist$gene_name,
                          keyType       = "SYMBOL",
                          OrgDb         = "org.Mm.eg.db",
                          ont           = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)
# visualize as dot plot
dotplot(bl4_enrichGO,showCategory = 10) 
