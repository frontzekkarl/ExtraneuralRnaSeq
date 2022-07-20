
# intersect terminal spleen and muscle changes with 3_vs_27 analysis from tabula muris senis
# import q<0.05 results, see above, for spleen and muscle
spleen273 <- read.delim("TabulaMurisSenisSpleen27vs3.csv",sep=",")
muscle273 <- read.delim("/home/ubuntu/Bulk_RNAseq/TabulaMurisSenisMuscle27vs3.csv",sep=",")
# separate up- and downregulated genes

spleen273_up <- spleen273[spleen273$log2FoldChange>=0.5,]
spleen273_down <- spleen273[spleen273$log2FoldChange<=-0.5,]
muscle273_up <- muscle273[muscle273$log2FoldChange>=0.5,]
muscle273_down <- muscle273[muscle273$log2FoldChange<=-0.5,]

sum(sp_term_up$gene_name %in% spleen273_up$gene_symbol) # 26
sum(sp_term_down$gene_name %in% spleen273_down$gene_symbol) # 1

sum(mu_term_up$gene_name %in% muscle273_up$gene_symbol) # 47
sum(mu_term_down$gene_name %in% muscle273_down$gene_symbol) # 96

## hypergeometric test for enrichment
# determine all genes expressed #
spleen273All <- read.delim("/home/ubuntu/Bulk_RNAseq/TabulaMurisSenisSpleen27vs3All.csv",sep=",")
muscle273All <- read.delim("/home/ubuntu/Bulk_RNAseq/TabulaMurisSenisMuscle27vs3All.csv",sep=",")
sp_termAll <- read.delim("/home/ubuntu/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/Spleen/result--RML6_term--over--NBH_term.txt")
mu_termAll <- read.delim("/home/ubuntu/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/MusclePlusBatch/result--RML6_term--over--NBH_term.txt")


spTermPresent <- sp_termAll[sp_termAll$isPresent==T,]
length(spTermPresent$gene_name) # 17360
# remove NA and logFoldChange=0 from spleen273All
sum(spleen273All$log2FoldChange!=0,na.rm=T) #23555

spleen273AllExpressed <- spleen273All[spleen273All$log2FoldChange!=0,] # select only expressed genes
sum(spleen273AllExpressed$gene_symbol %in% spTermPresent$gene_name) #13911
# total amount expressed spleen genes
17360+23555-13911 # 27004

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
## gmp not wokring in r studio, do calcs in cmd

# spleen maca vs. term
# total expressed genes, see above: 27004
# maca deg spleen 27vs3 w log2fc >= |0.5|: -> 283+15=298
length(spleen273_up$baseMean) #283
length(spleen273_down$baseMean) #15

# total deg spleen term: 692 up + 1753 down
# shared deg:
sum(spleen273_down$gene_symbol %in% sp_term_down$gene_name) #1
sum(spleen273_up$gene_symbol %in% sp_term_up$gene_name) #26

enrich_pvalue(27004,298,2445,27) # p=0.73
692+1753
# muscle

muTermPresent <- mu_termAll[mu_termAll$isPresent==T,]
length(muTermPresent$gene_name) # 14746, number of all expressed genes in muscle terminal
# determine all non-log2FoldChange==0 and non-NA in muscle273All
muscle273AllExpressed <- muscle273All[muscle273All$log2FoldChange!=0,]
# total muscle expressed genes
length(muscle273AllExpressed$gene_symbol) # 24061
# shared expressed genes
sum(muscle273AllExpressed$gene_symbol %in% muTermPresent$gene_name) #12647

# total amount expressed muscle genes
14746+24061-12647 # 26160

# deg muscle maca 27vs3
length(muscle273_down$baseMean) #461
length(muscle273_up$baseMean)#741

# deg muscle term: 612 up, 708 down=1320

# shared deg:
sum(muscle273_down$gene_symbol %in% mu_term_down$gene_name) #96
sum(muscle273_up$gene_symbol %in% mu_term_up$gene_name) #46 # one gene "Gcat" is twice in
# mu_term_up, was deleted -> 95
96+46
enrich_pvalue(26160,741,1320,141)  # 3.121034e-30

mu_term_up[mu_term_up$gene_name %in% muscle273_up$gene_symbol,]

# overlap 9vs1 months
spleen9vs1 <- read.delim("TabulaMurisSenisSpleen9vs1Sig.csv",sep=",")
muscle9vs1 <- read.delim("TabulaMurisSenisMuscle9vs1Sig.csv",sep=",")
#### Q<0.05 but LOG2 > 0.5 missing
#select up- and 3downregulated
spleen9vs1_up <- spleen9vs1[spleen9vs1$log2FoldChange>0,]
spleen9vs1_down <- spleen9vs1[spleen9vs1$log2FoldChange< 0,]

sum(sp_36_up$gene_name %in% spleen9vs1_up$gene_symbol) # 69 (310 deg sp_36UP,784 deg spleen9vs1UP)
sum(sp_36_down$gene_name %in% spleen9vs1_down$gene_symbol) #17 (55 sp36down, 226 spleen9vs1down)

### overlap prion deg w 27vs3 over time
## spleen
# 4wpi
sum(sp_4_down$gene_name %in% spleen273_down$gene_symbol) #0
# 16wpi
sum(sp_16_down$gene_name %in% spleen273_down$gene_symbol) #0
sum(sp_16_up$gene_name %in% spleen273_up$gene_symbol) #0

## muscle
# 4wpi
sum(mu_4_up$gene_name %in% muscle273_up$gene_symbol) #0
sum(mu_4_down$gene_name %in% muscle273_down$gene_symbol) #0


# 8wpi
sum(mu_8_up$gene_name %in% muscle273_up$gene_symbol) #0
sum(mu_8_down$gene_name %in% muscle273_down$gene_symbol) #1


#12wpi
sum(mu_12_up$gene_name %in% muscle273_up$gene_symbol) #1
sum(mu_12_down$gene_name %in% muscle273_down$gene_symbol) #2


# 14wpi
sum(mu_14_up$gene_name %in% muscle273_up$gene_symbol) #12
sum(mu_14_down$gene_name %in% muscle273_down$gene_symbol) #32

# 16wpi
sum(mu_16_up$gene_name %in% muscle273_up$gene_symbol) #0
sum(mu_16_down$gene_name %in% muscle273_down$gene_symbol) #0

# 18wpi
sum(mu_18_up$gene_name %in% muscle273_up$gene_symbol) #3
sum(mu_18_down$gene_name %in% muscle273_down$gene_symbol) #2

# 20wpi
sum(mu_20_up$gene_name %in% muscle273_up$gene_symbol) #1
sum(mu_20_down$gene_name %in% muscle273_down$gene_symbol) #13

# term
sum(mu_term_up$gene_name %in% muscle273_up$gene_symbol) #47
sum(mu_term_down$gene_name %in% muscle273_down$gene_symbol) #96 

## test mu log2>= |0.5|
muscle273_up <- muscle273[muscle273$log2FoldChange>=0.5,]
muscle273_down <- muscle273[muscle273$log2FoldChange<=-0.5,]

#### blood
bl273 <- read.delim("/home/ubuntu/Bulk_RNAseq/macabl273.csv",sep=",")
sum(bl273$gene_symbol %in%  bl_term$gene_name) #0

# muscle senescence
setwd("/home/ubuntu/Bulk_RNAseq/")
senes <- read.delim("unique_senescence_genes_murine.csv",sep=",")
muShared[muShared %in% senes$MGI.symbol]
mu_term[mu_term$gene_name == "Ets1",]
## ets1 is downregulated but this would most likely lead to protection against senescence
### 662      Ets1 0.0001938 0.004318    -0.5633
## DOI: 10.1126/sciadv.abf2017


senesHu <- read.delim("senescence_genes.csv",sep=";")

##GO Muscle term >>>> NOT WORKING <<<<<<<
# assemble list of shared genes
setwd("/home/ubuntu/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/")
muTermXprsd <- mu_termAll[mu_termAll$isPresent==T,]
muscle273 <- data.frame(results_list_MACA_bulk$Limb_Muscle$`27_vs_3`$ressig)
muShared <- muscle273_up$gene_symbol[muscle273_up$gene_symbol %in% mu_term_up$gene_name]
muSharedTemp <- muscle273_down$gene_symbol[muscle273_down$gene_symbol %in% mu_term_down$gene_name]
muShared <- append(muShared,muSharedTemp)
muTermXprsd <- mu_term$gene_name
write.csv(muShared,"/home/ubuntu/Bulk_RNAseq/muShared.csv",row.names = F)
muTermXprsd
library(org.Mm.eg.db)
muTermMaca <-           enrichGO(gene          = muShared,
                         universe      = muTermXprsd,
                         keyType       = "SYMBOL",
                         OrgDb         = "org.Mm.eg.db",
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
png(filename = 'mu',width=600 ,height = 300)
barplot(muTermMaca,showCategory = 5,col="p.adjust",font.size=12,border=TRUE,title='Blood, 4 wpi, upregulated \t [x = count]') 

# export shared muscle genes
muTermShared <- mu_term[mu_term$gene_name %in% muShared,]
mu273Shared <- muscle273[muscle273$gene_symbol %in% muShared,]
mu273SharedSorted <- mu273Shared[sort(mu273Shared$gene_symbol,index.return=T)$ix,]
muTermSharedSorted <- muTermShared[sort(muTermShared$gene_name,index.return=T)$ix,]
write.csv(muTermSharedSorted,"/home/ubuntu/Bulk_RNAseq/muTermSharedWith273.Sorted.csv")
write.csv(mu273SharedSorted,"/home/ubuntu/Bulk_RNAseq/mu273SharedWithTerm.Sorted.csv")

# export shared spleen genes
spTermShared <- sp_term[sp_term$gene_name %in% spShared,]
mu273Shared <- muscle273[muscle273$gene_symbol %in% muShared,]
mu273SharedSorted <- mu273Shared[sort(mu273Shared$gene_symbol,index.return=T)$ix,]
muTermSharedSorted <- muTermShared[sort(muTermShared$gene_name,index.return=T)$ix,]
