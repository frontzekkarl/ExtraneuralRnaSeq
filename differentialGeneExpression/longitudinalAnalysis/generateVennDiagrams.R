##
## create VennDiagrams for all organs and timepoints
## -------------------------------------------------

## 4 wpi
########
setwd("/home/ubuntu/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/")
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(patchwork)
library(UpSetR)

sp_4 <- read.delim("Spleen/result--RML6_4--over--NBH_4.txt")
bl_4 <- read.delim("BloodPlusBatch/result--RML6_4--over--NBH_4.txt")
mu_4 <- read.delim("MusclePlusBatch/result--RML6_4--over--NBH_4.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
brain_4 <- subset(brain_all,select=c(1,4:7))

# make lists with fdr<0.05 and |log2fc|>0.5
sp_4 <- subset(sp_4,(sp_4$fdr<0.05) & (abs(sp_4$log2.Ratio)>0.5),c(3,13,14,21))
bl_4 <- subset(bl_4,(bl_4$fdr<0.05) & (abs(bl_4$log2.Ratio)>0.5),c(3,13,14,21))
mu_4 <- subset(mu_4,(mu_4$fdr<0.05) & (abs(mu_4$log2.Ratio)>0.5),c(3,13,14,21))
br_4 <- subset(brain_4,(brain_4$Expressed_4wpi_main==TRUE) & (brain_4$FDR_4wpi_main<0.05) & (abs(brain_4$log2FC_4wpi_main)>0.5))
sp_4_up <- subset(sp_4,sp_4$log2.Ratio>0.5)
sp_4_down <- subset(sp_4,sp_4$log2.Ratio < 0)
br_4_up <- subset(br_4,br_4$log2FC_4wpi_main>0.5)
br_4_down <- subset(br_4,br_4$log2FC_4wpi_main<0)
bl_4_up <- subset(bl_4,bl_4$log2.Ratio>0.5)
bl_4_down <- subset(bl_4,bl_4$log2.Ratio<0)
mu_4_up <- subset(mu_4,mu_4$log2.Ratio>0.5)
mu_4_down <- subset(mu_4,mu_4$log2.Ratio<(-0.5))

library(VennDiagram)
library(RVenn)

v1 <- Venn(list(as.character(br_4_up$GeneID), 
                as.character(bl_4_up$gene_name), 
                as.character(sp_4_up$gene_name), 
                as.character(mu_4_up$gene_name)))
ggvenn(v1)
v1 <- venn.diagram(
  list(Brain=br_4_up$GeneID, Blood=bl_4_up$gene_name, Spleen=sp_4_up$gene_name, Muscle=mu_4_up$gene_name), 
  filename="Analysis/Venn/Venn_DEG_4wpi_up.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')

v2 <- venn.diagram(
  list(Brain=br_4_down$GeneID, Blood=bl_4_down$gene_name, Spleen=sp_4_down$gene_name, Muscle=mu_4_down$gene_name),
  filename="Analysis/Venn/Venn_DEG_4wpi_down.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')


## 8 wpi
########


sp_8 <- read.delim("Spleen/result--RML6_8--over--NBH_8.txt")
bl_8 <- read.delim("BloodPlusBatch/result--RML6_8--over--NBH_8.txt")
mu_8 <- read.delim("MusclePlusBatch/result--RML6_8--over--NBH_8.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
#colnames(brain_all)
brain_8 <- subset(brain_all,select=c(1,8:11))

# make lists with fdr<0.05 and |log2fc|>0.5
sp_8 <- subset(sp_8,(sp_8$fdr<0.05) & (abs(sp_8$log2.Ratio)>0.5),c(3,13,14,21))
bl_8 <- subset(bl_8,(bl_8$fdr<0.05) & (abs(bl_8$log2.Ratio)>0.5),c(3,13,14,21))
mu_8 <- subset(mu_8,(mu_8$fdr<0.05) & (abs(mu_8$log2.Ratio)>0.5),c(3,13,14,21))
br_8 <- subset(brain_8,(brain_8$Expressed_8wpi_main==TRUE) & (brain_8$FDR_8wpi_main<0.05) & (abs(brain_8$log2FC_8wpi_main)>0.5))
sp_8_up <- subset(sp_8,sp_8$log2.Ratio>0.5)
sp_8_down <- subset(sp_8,sp_8$log2.Ratio < 0)
br_8_up <- subset(br_8,br_8$log2FC_8wpi_main>0.5)
br_8_down <- subset(br_8,br_8$log2FC_8wpi_main<0)
bl_8_up <- subset(bl_8,bl_8$log2.Ratio>0.5)
bl_8_down <- subset(bl_8,bl_8$log2.Ratio<0)
mu_8_up <- subset(mu_8,mu_8$log2.Ratio>0.5)
mu_8_down <- subset(mu_8,mu_8$log2.Ratio<(-0.5))

v1 <- venn.diagram(
  list(Brain=br_8_up$GeneID, Blood=bl_8_up$gene_name, Spleen=sp_8_up$gene_name, Muscle=mu_8_up$gene_name), 
  filename="Analysis/Venn/Venn_DEG_8wpi_up.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')

v2 <- venn.diagram(
  list(Brain=br_8_down$GeneID, Blood=bl_8_down$gene_name, Spleen=sp_8_down$gene_name, Muscle=mu_8_down$gene_name),
  filename="Analysis/Venn/Venn_DEG_8wpi_down.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')


## 12 wpi
########


sp_12 <- read.delim("Spleen/result--RML6_12--over--NBH_12.txt")
bl_12 <- read.delim("BloodPlusBatch/result--RML6_12--over--NBH_12.txt")
mu_12 <- read.delim("MusclePlusBatch/result--RML6_12--over--NBH_12.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
colnames(brain_all)
brain_12 <- subset(brain_all,select=c(1,12:15))

# make lists with fdr<0.05 and |log2fc|>0.5
sp_12 <- subset(sp_12,(sp_12$fdr<0.05) & (abs(sp_12$log2.Ratio)>0.5),c(3,13,14,21))
bl_12 <- subset(bl_12,(bl_12$fdr<0.05) & (abs(bl_12$log2.Ratio)>0.5),c(3,13,14,21))
mu_12 <- subset(mu_12,(mu_12$fdr<0.05) & (abs(mu_12$log2.Ratio)>0.5),c(3,13,14,21))
br_12 <- subset(brain_12,(brain_12$Expressed_12wpi_main==TRUE) & (brain_12$FDR_12wpi_main<0.05) & (abs(brain_12$log2FC_12wpi_main)>0.5))
sp_12_up <- subset(sp_12,sp_12$log2.Ratio>0.5)
sp_12_down <- subset(sp_12,sp_12$log2.Ratio < 0)
br_12_up <- subset(br_12,br_12$log2FC_12wpi_main>0.5)
br_12_down <- subset(br_12,br_12$log2FC_12wpi_main<0)
bl_12_up <- subset(bl_12,bl_12$log2.Ratio>0.5)
bl_12_down <- subset(bl_12,bl_12$log2.Ratio<0)
mu_12_up <- subset(mu_12,mu_12$log2.Ratio>0.5)
mu_12_down <- subset(mu_12,mu_12$log2.Ratio<(-0.5))

br_12$GeneID[br_12$GeneID %in% mu_12$gene_name]
br_12[br_12$GeneID %in% br_12$GeneID[br_12$GeneID %in% mu_12$gene_name],]
mu_12[mu_12$gene_name %in% br_12$GeneID[br_12$GeneID %in% mu_12$gene_name],]


v1 <- venn.diagram(
  list(Brain=br_12_up$GeneID, Blood=bl_12_up$gene_name, Spleen=sp_12_up$gene_name, Muscle=mu_12_up$gene_name), 
  filename="Analysis/Venn/Venn_DEG_12wpi_up.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')

v2 <- venn.diagram(
  list(Brain=br_12_down$GeneID, Blood=bl_12_down$gene_name, Spleen=sp_12_down$gene_name, Muscle=mu_12_down$gene_name),
  filename="Analysis/Venn/Venn_DEG_12wpi_down.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')
    

## 14 wpi
#########
sp_14 <- subset(sp_14,(sp_14$fdr<0.05) & (abs(sp_14$log2.Ratio)>0.5),c(3,13,14,21))
bl_14 <- subset(bl_14,(bl_14$fdr<0.05) & (abs(bl_14$log2.Ratio)>0.5),c(3,13,14,21))
mu_14 <- subset(mu_14,(mu_14$fdr<0.05) & (abs(mu_14$log2.Ratio)>0.5),c(3,13,14,21))
br_14 <- subset(brain_14,(brain_14$Expressed_14wpi_main==TRUE) & (brain_14$FDR_14wpi_main<0.05) & (abs(brain_14$log2FC_14wpi_main)>0.5))
sp_14_up <- subset(sp_14,sp_14$log2.Ratio>0.5)
sp_14_down <- subset(sp_14,sp_14$log2.Ratio < 0)
br_14_up <- subset(br_14,br_14$log2FC_14wpi_main>0.5)
br_14_down <- subset(br_14,br_14$log2FC_14wpi_main<0)
bl_14_up <- subset(bl_14,bl_14$log2.Ratio>0.5)
bl_14_down <- subset(bl_14,bl_14$log2.Ratio<0)
mu_14_up <- subset(mu_14,mu_14$log2.Ratio>0.5)
mu_14_down <- subset(mu_14,mu_14$log2.Ratio<(-0.5))

br_14$GeneID[br_14$GeneID %in% mu_14$gene_name]
br_14[br_14$GeneID %in% br_14$GeneID[br_14$GeneID %in% mu_14$gene_name],]
mu_14[mu_14$gene_name %in% br_14$GeneID[br_14$GeneID %in% mu_14$gene_name],]

v1 <- venn.diagram(
  list(Brain=br_14_up$GeneID, Blood=bl_14_up$gene_name, Spleen=sp_14_up$gene_name, Muscle=mu_14_up$gene_name), 
  filename="Analysis/Venn/Venn_DEG_14wpi_up.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')

v2 <- venn.diagram(
  list(Brain=br_14_down$GeneID, Blood=bl_14_down$gene_name, Spleen=sp_14_down$gene_name, Muscle=mu_14_down$gene_name),
  filename="Analysis/Venn/Venn_DEG_14wpi_down.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')



## 16 wpi
########


sp_16 <- read.delim("Spleen/result--RML6_16--over--NBH_16.txt")
bl_16 <- read.delim("BloodPlusBatch/result--RML6_16--over--NBH_16.txt")
mu_16 <- read.delim("MusclePlusBatch/result--RML6_16--over--NBH_16.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
colnames(brain_all)
brain_16 <- subset(brain_all,select=c(1,20:24))

# make lists with fdr<0.05 and |log2fc|>0.5
sp_16 <- subset(sp_16,(sp_16$fdr<0.05) & (abs(sp_16$log2.Ratio)>0.5),c(3,13,14,21))
bl_16 <- subset(bl_16,(bl_16$fdr<0.05) & (abs(bl_16$log2.Ratio)>0.5),c(3,13,14,21))
mu_16 <- subset(mu_16,(mu_16$fdr<0.05) & (abs(mu_16$log2.Ratio)>0.5),c(3,13,14,21))
br_16 <- subset(brain_16,(brain_16$Expressed_16wpi_main==TRUE) & (brain_16$FDR_16wpi_main<0.05) & (abs(brain_16$log2FC_16wpi_main)>0.5))
sp_16_up <- subset(sp_16,sp_16$log2.Ratio>0.5)
sp_16_down <- subset(sp_16,sp_16$log2.Ratio < 0)
br_16_up <- subset(br_16,br_16$log2FC_16wpi_main>0.5)
br_16_down <- subset(br_16,br_16$log2FC_16wpi_main<0)
bl_16_up <- subset(bl_16,bl_16$log2.Ratio>0.5)
bl_16_down <- subset(bl_16,bl_16$log2.Ratio<0)
mu_16_up <- subset(mu_16,mu_16$log2.Ratio>0.5)
mu_16_down <- subset(mu_16,mu_16$log2.Ratio<(-0.5))

v1 <- venn.diagram(
  list(Brain=br_16_up$GeneID, Blood=bl_16_up$gene_name, Spleen=sp_16_up$gene_name, Muscle=mu_16_up$gene_name), 
  filename="Analysis/Venn/Venn_DEG_16wpi_up.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')

v2 <- venn.diagram(
  list(Brain=br_16_down$GeneID, Blood=bl_16_down$gene_name, Spleen=sp_16_down$gene_name, Muscle=mu_16_down$gene_name),
  filename="Analysis/Venn/Venn_DEG_16wpi_down.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')


## 18 wpi
########


sp_18 <- read.delim("Spleen/result--RML6_18--over--NBH_18.txt")
bl_18 <- read.delim("BloodPlusBatch/result--RML6_18--over--NBH_18.txt")
mu_18 <- read.delim("MusclePlusBatch/result--RML6_18--over--NBH_18.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
colnames(brain_all)
brain_18 <- subset(brain_all,select=c(1,24:27))
head(brain_18)
# make lists with fdr<0.05 and |log2fc|>0.5
sp_18 <- subset(sp_18,(sp_18$fdr<0.05) & (abs(sp_18$log2.Ratio)>0.5),c(3,13,14,21))
bl_18 <- subset(bl_18,(bl_18$fdr<0.05) & (abs(bl_18$log2.Ratio)>0.5),c(3,13,14,21))
mu_18 <- subset(mu_18,(mu_18$fdr<0.05) & (abs(mu_18$log2.Ratio)>0.5),c(3,13,14,21))
br_18 <- subset(brain_18,(brain_18$Expressed_18wpi_main==TRUE) & (brain_18$FDR_18wpi_main<0.05) & (abs(brain_18$log2FC_18wpi_main)>0.5))
sp_18_up <- subset(sp_18,sp_18$log2.Ratio>0.5)
sp_18_down <- subset(sp_18,sp_18$log2.Ratio < 0)
br_18_up <- subset(br_18,br_18$log2FC_18wpi_main>0.5)
br_18_down <- subset(br_18,br_18$log2FC_18wpi_main<0)
bl_18_up <- subset(bl_18,bl_18$log2.Ratio>0.5)
bl_18_down <- subset(bl_18,bl_18$log2.Ratio<0)
mu_18_up <- subset(mu_18,mu_18$log2.Ratio>0.5)
mu_18_down <- subset(mu_18,mu_18$log2.Ratio<(-0.5))

v1 <- venn.diagram(
  list(Brain=br_18_up$GeneID, Blood=bl_18_up$gene_name, Spleen=sp_18_up$gene_name, Muscle=mu_18_up$gene_name), 
  filename="Analysis/Venn/Venn_DEG_18wpi_up.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')

v2 <- venn.diagram(
  list(Brain=br_18_down$GeneID, Blood=bl_18_down$gene_name, Spleen=sp_18_down$gene_name, Muscle=mu_18_down$gene_name),
  filename="Analysis/Venn/Venn_DEG_18wpi_down.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')



## 20 wpi
########

sp_20 <- read.delim("Spleen/result--RML6_20--over--NBH_20.txt")
bl_20 <- read.delim("BloodPlusBatch/result--RML6_20--over--NBH_20.txt")
mu_20 <- read.delim("MusclePlusBatch/result--RML6_20--over--NBH_20.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
colnames(brain_all)
brain_20 <- subset(brain_all,select=c(1,28:31))
head(brain_20)
# make lists with fdr<0.05 and |log2fc|>0.5
sp_20 <- subset(sp_20,(sp_20$fdr<0.05) & (abs(sp_20$log2.Ratio)>0.5),c(3,13,14,21))
bl_20 <- subset(bl_20,(bl_20$fdr<0.05) & (abs(bl_20$log2.Ratio)>0.5),c(3,13,14,21))
mu_20 <- subset(mu_20,(mu_20$fdr<0.05) & (abs(mu_20$log2.Ratio)>0.5),c(3,13,14,21))
br_20 <- subset(brain_20,(brain_20$Expressed_20wpi_main==TRUE) & (brain_20$FDR_20wpi_main<0.05) & (abs(brain_20$log2FC_20wpi_main)>0.5))
sp_20_up <- subset(sp_20,sp_20$log2.Ratio>0.5)
sp_20_down <- subset(sp_20,sp_20$log2.Ratio < 0)
br_20_up <- subset(br_20,br_20$log2FC_20wpi_main>0.5)
br_20_down <- subset(br_20,br_20$log2FC_20wpi_main<0)
bl_20_up <- subset(bl_20,bl_20$log2.Ratio>0.5)
bl_20_down <- subset(bl_20,bl_20$log2.Ratio<0)
mu_20_up <- subset(mu_20,mu_20$log2.Ratio>0.5)
mu_20_down <- subset(mu_20,mu_20$log2.Ratio<(-0.5))

v1 <- venn.diagram(
  list(Brain=br_20_up$GeneID, Blood=bl_20_up$gene_name, Spleen=sp_20_up$gene_name, Muscle=mu_20_up$gene_name), 
  filename="Analysis/Venn/Venn_DEG_20wpi_up.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')

v2 <- venn.diagram(
  list(Brain=br_20_down$GeneID, Blood=bl_20_down$gene_name, Spleen=sp_20_down$gene_name, Muscle=mu_20_down$gene_name),
  filename="Analysis/Venn/Venn_DEG_20wpi_down.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')



## terminal disease
########


sp_term <- read.delim("Spleen/result--RML6_term--over--NBH_term.txt")
bl_term <- read.delim("BloodPlusBatch/result--RML6_term--over--NBH_term.txt")
mu_term <- read.delim("MusclePlusBatch/result--RML6_term--over--NBH_term.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
brain_term <- subset(brain_all,select=c(1,32:35))

# make lists with fdr<0.05 and |log2fc|>0.5
sp_term <- subset(sp_term,(sp_term$fdr<0.05) & (abs(sp_term$log2.Ratio)>0.5),c(3,13,14,21))
bl_term <- subset(bl_term,(bl_term$fdr<0.05) & (abs(bl_term$log2.Ratio)>0.5),c(3,13,14,21))
mu_term <- subset(mu_term,(mu_term$fdr<0.05) & (abs(mu_term$log2.Ratio)>0.5),c(3,13,14,21))
br_term <- subset(brain_term,(brain_term$Expressed_term_main==TRUE) & (brain_term$FDR_term_main<0.05) & (abs(brain_term$log2FC_term_main)>0.5))
sp_term_up <- subset(sp_term,sp_term$log2.Ratio>0.5)
sp_term_down <- subset(sp_term,sp_term$log2.Ratio < 0)
br_term_up <- subset(br_term,br_term$log2FC_term_main>0.5)
br_term_down <- subset(br_term,br_term$log2FC_term_main<0)
bl_term_up <- subset(bl_term,bl_term$log2.Ratio>0.5)
bl_term_down <- subset(bl_term,bl_term$log2.Ratio<0)
mu_term_up <- subset(mu_term,mu_term$log2.Ratio>0.5)
mu_term_down <- subset(mu_term,mu_term$log2.Ratio<(-0.5))


v1 <- venn.diagram(
  list(Brain=br_term_up$GeneID, Blood=bl_term_up$gene_name, Spleen=sp_term_up$gene_name, Muscle=mu_term_up$gene_name), 
  filename="Analysis/Venn_DEG_terminal_up.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')

v2 <- venn.diagram(
  list(Brain=br_term_down$GeneID, Blood=bl_term_down$gene_name, Spleen=sp_term_down$gene_name, Muscle=mu_term_down$gene_name),
  filename="Analysis/Venn_DEG_terminal_down.png", 
  height=5000,
  width=5000,
  cex=4,
  cat.cex=4,
  fill=rainbow(4),
  main.fontfamily='arial',
  fontfamily='arial')
