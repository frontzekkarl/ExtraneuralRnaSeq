
# intersect terminal spleen and muscle changes with 3_vs_27 analysis from tabula muris senis, have to be download prior to analysis from AWS:
# https://registry.opendata.aws/tabula-muris-senis/
# import q<0.05 results, see above, for spleen and muscle
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
# intersect terminal spleen and muscle changes with 3_vs_27 analysis from tabula muris senis
# import q<0.05 results, see above, for spleen and muscle
spleen273 <- read.delim("TabulaMurisSenisSpleen27vs3.csv",sep=",")
muscle273 <- read.delim("/home/ubuntu/Bulk_RNAseq/TabulaMurisSenisMuscle27vs3.csv",sep=",")

# separate up- and downregulated genes

spleen273_up <- spleen273[spleen273$log2FoldChange>=0.5,]
spleen273_down <- spleen273[spleen273$log2FoldChange<=-0.5,]
muscle273_up <- muscle273[muscle273$log2FoldChange>=0.5,]
muscle273_down <- muscle273[muscle273$log2FoldChange<=-0.5,]

# subset terminal degs, choose fdr<0.05, |log2fc|>0.5


sp_term <- read.delim("/home/ubuntu/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/Spleen/result--RML6_term--over--NBH_term.txt")
mu_term <- read.delim("/home/ubuntu/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/MusclePlusBatch/result--RML6_term--over--NBH_term.txt")
sp_term <- subset(sp_term,(sp_term$fdr<0.05) & (abs(sp_term$log2.Ratio)>0.5),c(3,13,14,21))
mu_term <- subset(mu_term,(mu_term$fdr<0.05) & (abs(mu_term$log2.Ratio)>0.5),c(3,13,14,21))
sp_term_up <- subset(sp_term,sp_term$log2.Ratio>0.5)
sp_term_down <- subset(sp_term,sp_term$log2.Ratio < 0)
mu_term_up <- subset(mu_term,mu_term$log2.Ratio>0.5)
mu_term_down <- subset(mu_term,mu_term$log2.Ratio<(-0.5))


# count overlap between groups
sum(sp_term_up$gene_name %in% spleen273_up$gene_symbol) # 26
sum(sp_term_down$gene_name %in% spleen273_down$gene_symbol) # 1

sum(mu_term_up$gene_name %in% muscle273_up$gene_symbol) # 39
sum(mu_term_down$gene_name %in% muscle273_down$gene_symbol) # 88

## hypergeometric test for enrichment
# determine all genes expressesd #
spleen273All <- read.delim("/home/ubuntu/Bulk_RNAseq/TabulaMurisSenisSpleen27vs3All.csv",sep=",")
muscle273All <- read.delim("/home/ubuntu/Bulk_RNAseq/TabulaMurisSenisMuscle27vs3All.csv",sep=",")


##### calculate overlap for splenic genes
# reminder of overlapping genes
sum(sp_term_up$gene_name %in% spleen273_up$gene_symbol) # 26
sum(sp_term_down$gene_name %in% spleen273_down$gene_symbol) # 1
#total overlapping genes =27

# determine expressed genes in terminal spleen
# read sp_term2 for field "isPresent"
sp_term2 <- read.delim("/home/ubuntu/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/Spleen/result--RML6_term--over--NBH_term.txt")
sp_term_present <- subset(sp_term2,(sp_term2$isPresent==TRUE))

# calculate overlap of expressed genes with tabula rasa, used as background space
sum(sp_term_present$gene_name %in% spleen273All$gene_symbol) #14116

##### calculate total DEG per group

length(spleen273_up$baseMean) #283
length(spleen273_down$baseMean) #15
# total in spleen tabula senis = 298

length(sp_term_up$gene_name) # 692
length(sp_term_down$gene_name) # 1753
#total in spleen prion = 2445

enrich_pvalue(14116,311,2445,29) # P 0.99


### same for muscle

# shared deg between prion and tabula muris senis
sum(mu_term_up$gene_name %in% muscle273_up$gene_symbol) # 39
sum(mu_term_down$gene_name %in% muscle273_down$gene_symbol) # 88

# total deg tabula muris: 208+461=741
length(muscle273_up$baseMean) # 280
length(muscle273_down$baseMean) # 461

# same as above, read mu_term2 for field #isPresent
mu_term2 <- read.delim("/home/ubuntu/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/MusclePlusBatch/result--RML6_term--over--NBH_term.txt")
mu_term_present <- subset(mu_term2,mu_term2$isPresent==TRUE) # subset present genes

length(mu_term_up$gene_name) # 612
length(mu_term_down$gene_name) # 708
# total amount deg muscle 1320

# generate backgroudn expression set
sum(mu_term_present$gene_name %in% muscle273All$gene_symbol) #12774


enrich_pvalue(12774,669,1320,39+88) # P 0.0012

#####
#####shared genes Pmid29382830 vs blood terminal
# read bl_term Fdr and log2 fold adjusted from shared_genes_timecourse.R

setwd("~/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/BloodPlusBatch")
bl_term_raw <- read.delim("result--RML6_term--over--NBH_term.txt")

bl_term_raw <- bl_term_raw[bl_term_raw$isPresent==TRUE,]

setwd("~/Bulk_RNAseq/Pmid29382830")
bl_aged<- read.delim("blood_9_30_edgeR_extended.csv",header=F)
# count total expressed genes shared by Pmid29382830 and blood terminal

colnames(bl_aged) <- c("geneid",'log2fc','log2cpm','pval','fdr','ensg','geneid2','desc','func','annot')
bl_aged <- subset(bl_aged,abs(bl_aged$log2fc)>0)
sum(bl_term_raw$gene_name %in% bl_aged$geneid) #4475

# select deg at log2fc > 0.5 and fdr <0.05
bl_aged <- subset(bl_aged,(abs(bl_aged$log2fc)>0.5) & (bl_aged$fdr<0.05), select=c("geneid",'log2fc','pval','fdr'))
bl_term <- subset(bl_term_raw,(bl_term_raw$fdr<0.05) & (abs(bl_term_raw$log2.Ratio)>0.5),c(3,13,14,21))

bl_term_up <- subset(bl_term,bl_term$log2.Ratio>0.5)
bl_term_down <- subset(bl_term,bl_term$log2.Ratio<0)

length(bl_term_up$gene_name)#5
length(bl_term_down$gene_name)#11
# total deg blood prion 16

bl_aged_up <- subset(bl_aged,bl_aged$log2fc>0.5)
bl_aged_down <- subset(bl_aged,bl_aged$log2fc<0)
length(bl_aged_up$geneid) #362
length(bl_aged_down$geneid) #198
# total deg blood aged mice 560

# test for enrichment
# gene expression space 4475
# total deg aged 560
# total deg term 16
# shared 1

enrich_pvalue(4475,560,16,1) # 0.90
