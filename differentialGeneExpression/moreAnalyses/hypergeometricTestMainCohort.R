## test enrichment of DEGs between main cohort organs and timepoints

##### read DEGs to determine total amount of DEGs as well as shared DEGs
sp_4 <- read.delim("Spleen/result--RML6_4--over--NBH_4.txt")
bl_4 <- read.delim("BloodPlusBatch/result--RML6_4--over--NBH_4.txt")
mu_4 <- read.delim("MusclePlusBatch/result--RML6_4--over--NBH_4.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
brain_4 <- subset(brain_all,select=c(1,4:7))

sp_8 <- read.delim("Spleen/result--RML6_8--over--NBH_8.txt")
bl_8 <- read.delim("BloodPlusBatch/result--RML6_8--over--NBH_8.txt")
mu_8 <- read.delim("MusclePlusBatch/result--RML6_8--over--NBH_8.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
#colnames(brain_all)
brain_8 <- subset(brain_all,select=c(1,8:11))

mu_12 <- read.delim("MusclePlusBatch/result--RML6_12--over--NBH_12.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
colnames(brain_all)
brain_12 <- subset(brain_all,select=c(1,12:15))

sp_14 <- read.delim("Spleen/result--RML6_14--over--NBH_14.txt")
bl_14 <- read.delim("BloodPlusBatch/result--RML6_14--over--NBH_14.txt")
mu_14 <- read.delim("MusclePlusBatch/result--RML6_14--over--NBH_14.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
colnames(brain_all)
brain_14 <- subset(brain_all,select=c(1,16:19))


sp_18 <- read.delim("Spleen/result--RML6_18--over--NBH_18.txt")
bl_18 <- read.delim("BloodPlusBatch/result--RML6_18--over--NBH_18.txt")
mu_18 <- read.delim("MusclePlusBatch/result--RML6_18--over--NBH_18.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
colnames(brain_all)
brain_18 <- subset(brain_all,select=c(1,24:27))

sp_20 <- read.delim("Spleen/result--RML6_20--over--NBH_20.txt")
bl_20 <- read.delim("BloodPlusBatch/result--RML6_20--over--NBH_20.txt")
mu_20 <- read.delim("MusclePlusBatch/result--RML6_20--over--NBH_20.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
colnames(brain_all)
brain_20 <- subset(brain_all,select=c(1,28:31))

sp_term <- read.delim("Spleen/result--RML6_term--over--NBH_term.txt")
bl_term <- read.delim("BloodPlusBatch/result--RML6_term--over--NBH_term.txt")
mu_term <- read.delim("MusclePlusBatch/result--RML6_term--over--NBH_term.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
brain_term <- subset(brain_all,select=c(1,32:35))

####### get list of all expressed genes in overlaps
### 4 wks
group1=bl_4
group2=mu_4
g1_expressed=sum(group1$isPresent == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$isPresent == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$gene_name %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared
brain_8

### 8 wks
group1=brain_8
group2=mu_8
g1_expressed=sum(group1$Expressed_8wpi_main == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$Expressed_8wpi_main == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$GeneID %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared



# 12 wpi muscle + brain
group1=brain_12
group2=mu_12
g1_expressed=sum(group1$Expressed_12wpi_main == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$Expressed_12wpi_main == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$GeneID %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared

# 14 wpi muscle + brain
group1=brain_14
group2=mu_14
g1_expressed=sum(group1$Expressed_14wpi_main == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$Expressed_14wpi_main == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$GeneID %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared

# 14 wpi muscle + blood
group1=bl_14
group2=mu_14
g1_expressed=sum(group1$isPresent == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$isPresent == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$gene_name %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared

# 18 wpi muscle + brain
group1=brain_18
group2=mu_18
g1_expressed=sum(group1$Expressed_18wpi_main == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$Expressed_18wpi_main == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$GeneID %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared

# 18 wpi muscle + blood
group1=bl_18
group2=mu_18
g1_expressed=sum(group1$isPresent == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$isPresent == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$gene_name %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared

# 18 wpi brain + blood
group1=brain_18
group2=bl_18
g1_expressed=sum(group1$Expressed_18wpi_main == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$Expressed_18wpi_main == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$GeneID %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared

# 20 wpi brain + blood
group1=brain_20
group2=bl_20
g1_expressed=sum(group1$Expressed_20wpi_main == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$Expressed_20wpi_main == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$GeneID %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared

# 20 wpi muscle + blood
group1=bl_20
group2=mu_20
g1_expressed=sum(group1$isPresent == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$isPresent == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$gene_name %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared

# term brain + blood
group1=brain_term
group2=mu_term
g1_expressed=sum(group1$Expressed_term_main == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$Expressed_term_main == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$GeneID %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared

# term brain + spleen
group1=brain_term
group2=sp_term
g1_expressed=sum(group1$Expressed_term_main == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$Expressed_term_main == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$GeneID %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared


# term muscle + spleen
group1=sp_term
group2=mu_term
g1_expressed=sum(group1$isPresent == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$isPresent == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$gene_name %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared


# term blood + spleen
group1=sp_term
group2=bl_term
g1_expressed=sum(group1$isPresent == TRUE)
g2_expressed=sum(group2$isPresent == TRUE)
group1_xg=subset(group1, group1$isPresent == TRUE)
group2_xg=subset(group2, group2$isPresent == TRUE)
g12_shared=sum(group1_xg$gene_name %in% group2_xg$gene_name) # 10876
# number of all genes -shared
g1_expressed+g2_expressed-g12_shared


# hypergeometric test function
# N : all expressed genes between two experimental groups
# A : DEGs from group1
# B : DEGs from group2
# k : shared DEGs

library(gmp) 
enrich_pvalue <- function(N, A, B, k)
{
  m <- A + k
  n <- B + k
  i <- k:min(m,n)
  
  as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
}


# 4 wpi blood + muscle
enrich_pvalue(16860,319,26,1) # p=0.40
# 8 wpi muscle + brain
enrich_pvalue(15556,95+11,101+712,9) # p=0.15
# 12 wpi muscle + brain
enrich_pvalue(15737,36+19,143+80,5) # p=0.0017
# 14 wpi muscle + brain
enrich_pvalue(14375,11+77,189+287,5) # p=0.20
# 14 wpi muscle + blood
enrich_pvalue(16842,189+287,100,6) # p=0.084
# 18 wpi muscle + brain
enrich_pvalue(14384,211+80,338+24,17) # p=0.0033
# 18 wpi muscle + blood
enrich_pvalue(16468,138+63,211+80,1) # p=0.97
# 18 wpi brain + blood
enrich_pvalue(15675,211+80,338+24,2) # p=0.99
# 20 wpi muscle + brain 
enrich_pvalue(15814,120+77,467+146,8) # p=0.56
# 20 wpi muscle + blood
enrich_pvalue(17073,120+77,8,1) # p = 0.09
# terminal muscle + brain
enrich_pvalue(15974,1815+1013,612+708,56+54) # p = 1
# terminal spleen + brain
enrich_pvalue(18300,1815+1013,692+1753,74+39) # p = 1
# terminal spleen + muscle
enrich_pvalue(18524,612+708,692+1753,53+42) # p = 1
# terminal spleen + blood
enrich_pvalue(17688,16,692+1753,9) # p = 0.0047

### test for enriched genes

# muscle brain wpi8
mu_8 <- subset(mu_8,(mu_8$fdr<0.05) & (abs(mu_8$log2.Ratio)>0.5),c(3,13,14,21))
br_8 <- subset(brain_8,(brain_8$Expressed_8wpi_main==TRUE) & (brain_8$FDR_8wpi_main<0.05) & (abs(brain_8$log2FC_8wpi_main)>0.5))

mubr_8=mu_8[mu_8$gene_name %in% br_8$GeneID,]

# muscle brain wpi12
mu_12 <- subset(mu_12,(mu_12$fdr<0.05) & (abs(mu_12$log2.Ratio)>0.5),c(3,13,14,21))
br_12 <- subset(brain_12,(brain_12$Expressed_12wpi_main==TRUE) & (brain_12$FDR_12wpi_main<0.05) & (abs(brain_12$log2FC_12wpi_main)>0.5))

mubr_12=mu_12[mu_12$gene_name %in% br_12$GeneID,]

# muscle brain wpi14
mu_14 <- subset(mu_14,(mu_14$fdr<0.05) & (abs(mu_14$log2.Ratio)>0.5),c(3,13,14,21))
br_14 <- subset(brain_14,(brain_14$Expressed_14wpi_main==TRUE) & (brain_14$FDR_14wpi_main<0.05) & (abs(brain_14$log2FC_14wpi_main)>0.5))

mubr_14=mu_14[mu_14$gene_name %in% br_14$GeneID,]

mubr_12[mubr_12$gene_name %in% mubr_14$gene_name,]

# muscle brain wpi18
mu_18 <- subset(mu_18,(mu_18$fdr<0.05) & (abs(mu_18$log2.Ratio)>0.5),c(3,13,14,21))
br_18 <- subset(brain_18,(brain_18$Expressed_18wpi_main==TRUE) & (brain_18$FDR_18wpi_main<0.05) & (abs(brain_18$log2FC_18wpi_main)>0.5))

mubr_18=mu_18[mu_18$gene_name %in% br_18$GeneID,]
mubr_18$gene_name %in% mubr_14$gene_name

# muscle blood wpi14
bl_14 <- subset(bl_14,(bl_14$fdr<0.05) & (abs(bl_14$log2.Ratio)>0.5),c(3,13,14,21))
mubl_14=mu_14[mu_14$gene_name %in% bl_14$gene_name,]
mubl_14[mubl_14$gene_name %in% mubr_14$gene_name,]
mubr_14


# muscle blood wpi20
bl_20 <- subset(bl_20,(bl_20$fdr<0.05) & (abs(bl_20$log2.Ratio)>0.5),c(3,13,14,21))
mu_20 <- subset(mu_20,(mu_20$fdr<0.05) & (abs(mu_20$log2.Ratio)>0.5),c(3,13,14,21))

mubl_20=mu_20[mu_20$gene_name %in% bl_20$gene_name,]
mubl_20[mubl_20$gene_name %in% mubr_18$gene_name,]

# spleen blood term
sp_term <- subset(sp_term,(sp_term$fdr<0.05) & (abs(sp_term$log2.Ratio)>0.5),c(3,13,14,21))
bl_term <- subset(bl_term,(bl_term$fdr<0.05) & (abs(bl_term$log2.Ratio)>0.5),c(3,13,14,21))

blsp_term=sp_term[sp_term$gene_name %in% bl_term$gene_name,]
blsp_term$gene_name %in% mubr_8$gene_name


