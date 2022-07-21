# hierarchical clustering of muscle, main cohort
# read mu_4 - mu_term according to "analysisMuscleMainCohort.R" from /differentialGeneExpression

library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(pheatmap)

# sequential full join of log2 ratios of all degs at log2ratio > 0.5, fdr < 0.05
##########

# make gene_name > rowname

rownames(mu_4) <- mu_4$gene_name
rownames(mu_8) <- mu_8$gene_name

mu_12 <- mu_12[!duplicated(mu_12$gene_name),] # delete duplicates
rownames(mu_12) <- mu_12$gene_name

length(mu_14$gene_name) == length(unique(mu_14$gene_name)) # false
mu_14 <- mu_14[!duplicated(mu_14$gene_name),] # delete duplicates
rownames(mu_14) <- mu_14$gene_name

length(mu_16$gene_name) == length(unique(mu_16$gene_name)) # true
rownames(mu_16) <- mu_16$gene_name

length(mu_18$gene_name) == length(unique(mu_18$gene_name)) # false
mu_18 <- mu_18[!duplicated(mu_18$gene_name),] # delete duplicates
rownames(mu_18) <- mu_18$gene_name

length(mu_20$gene_name) == length(unique(mu_20$gene_name)) # true
rownames(mu_20) <- mu_20$gene_name

length(mu_term$gene_name) == length(unique(mu_term$gene_name)) # false
mu_term <- mu_term[!duplicated(mu_term$gene_name),] # delete duplicates
rownames(mu_term) <- mu_term$gene_name

mumu <- merge(mu_4,mu_8,by=0,all=T)
rownames(mumu) <- mumu$Row.names
mumu <- subset(mumu,select=c('log2.Ratio.x','log2.Ratio.y'))
colnames(mumu)=c('4wpi','8wpi')

mumu <- merge(mumu,mu_12,by=0,all=T)
rownames(mumu) <- mumu$Row.names
mumu <- subset(mumu,select=c('4wpi','8wpi','log2.Ratio'))
colnames(mumu)=c('4wpi','8wpi','12wpi')

mumu <- merge(mumu,mu_14,by=0,all=T)
rownames(mumu) <- mumu$Row.names
mumu <- subset(mumu,select=c('4wpi','8wpi','12wpi','log2.Ratio'))
colnames(mumu)=c('4wpi','8wpi','12wpi','14wpi')


mumu <- merge(mumu,mu_16,by=0,all=T)
rownames(mumu) <- mumu$Row.names
mumu <- subset(mumu,select=c('4wpi','8wpi','12wpi','14wpi','log2.Ratio'))
colnames(mumu)=c('4wpi','8wpi','12wpi','14wpi','16wpi')

mumu <- merge(mumu,mu_18,by=0,all=T)
rownames(mumu) <- mumu$Row.names
mumu <- subset(mumu,select=c('4wpi','8wpi','12wpi','14wpi','16wpi','log2.Ratio'))
colnames(mumu)=c('4wpi','8wpi','12wpi','14wpi','16wpi','18wpi')

mumu <- merge(mumu,mu_20,by=0,all=T)
rownames(mumu) <- mumu$Row.names
mumu <- subset(mumu,select=c('4wpi','8wpi','12wpi','14wpi','16wpi','18wpi','log2.Ratio'))
colnames(mumu)=c('4wpi','8wpi','12wpi','14wpi','16wpi','18wpi','20wpi')


mumu <- merge(mumu,mu_term,by=0,all=T)
rownames(mumu) <- mumu$Row.names
mumu <- subset(mumu,select=c('4wpi','8wpi','12wpi','14wpi','16wpi','18wpi','20wpi','log2.Ratio'))
colnames(mumu)=c('4wpi','8wpi','12wpi','14wpi','16wpi','18wpi','20wpi','term')

### plot heatmap


mumu <- data.matrix(mumu) # create data matrix for heatmap
mumu[is.na(mumu)] <- 0 # change na -> 0

pheatmap(mumu,
         cluster_cols=F,
         show_rownames = F,
         colorRampPalette(c("#0512B7", "white", "#E83DFC"))(75)
)

