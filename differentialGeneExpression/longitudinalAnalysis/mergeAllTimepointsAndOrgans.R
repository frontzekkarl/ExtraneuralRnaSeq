# read all files
setwd("/home/ubuntu/Bulk_RNAseq/Comparisons_NewPipeline/p3506_PeripheralSamples/")
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

sp_8 <- read.delim("Spleen/result--RML6_8--over--NBH_8.txt")
bl_8 <- read.delim("BloodPlusBatch/result--RML6_8--over--NBH_8.txt")
mu_8 <- read.delim("MusclePlusBatch/result--RML6_8--over--NBH_8.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
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




sp_12 <- read.delim("Spleen/result--RML6_12--over--NBH_12.txt")
bl_12 <- read.delim("BloodPlusBatch/result--RML6_12--over--NBH_12.txt")
mu_12 <- read.delim("MusclePlusBatch/result--RML6_12--over--NBH_12.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
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

sp_14 <- read.delim("Spleen/result--RML6_14--over--NBH_14.txt")
bl_14 <- read.delim("BloodPlusBatch/result--RML6_14--over--NBH_14.txt")
mu_14 <- read.delim("MusclePlusBatch/result--RML6_14--over--NBH_14.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
brain_14 <- subset(brain_all,select=c(1,16:19))

# make lists with fdr<0.05 and |log2fc|>0.5
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


sp_16 <- read.delim("Spleen/result--RML6_16--over--NBH_16.txt")
bl_16 <- read.delim("BloodPlusBatch/result--RML6_16--over--NBH_16.txt")
mu_16 <- read.delim("MusclePlusBatch/result--RML6_16--over--NBH_16.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
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

sp_18 <- read.delim("Spleen/result--RML6_18--over--NBH_18.txt")
bl_18 <- read.delim("BloodPlusBatch/result--RML6_18--over--NBH_18.txt")
mu_18 <- read.delim("MusclePlusBatch/result--RML6_18--over--NBH_18.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
brain_18 <- subset(brain_all,select=c(1,24:27))

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

sp_20 <- read.delim("Spleen/result--RML6_20--over--NBH_20.txt")
bl_20 <- read.delim("BloodPlusBatch/result--RML6_20--over--NBH_20.txt")
mu_20 <- read.delim("MusclePlusBatch/result--RML6_20--over--NBH_20.txt")
brain_all <- read.delim("Analysis/brain_plospathogens.txt")
brain_20 <- subset(brain_all,select=c(1,28:31))

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

# assign timepoints and organ names to tables prior to merging
sp_4_down$wpi <- 4
sp_4_down$organ <- "spleen"
sp_16_up$wpi <- 16
sp_16_up$organ <- "spleen"
sp_16_down$wpi <- 16
sp_16_down$organ <- "spleen"
sp_term_down$wpi <- 22
sp_term_down$organ <- "spleen"
sp_term_up$wpi <- 22
sp_term_up$organ <- "spleen"

bl_4_up$wpi <- 4
bl_4_down$wpi <- 4
bl_14_down$wpi <- 14
bl_16_down$wpi <- 16
bl_16_up$wpi <- 16
bl_18_up$wpi <- 18
bl_18_down$wpi <- 18
bl_20_down$wpi <- 20
bl_20_up$wpi <- 20
bl_term_up$wpi <- 22
bl_term_down$wpi <- 22
bl_4_up$organ <- "blood"
bl_4_down$organ <- "blood"
bl_14_down$organ <- "blood"
bl_16_down$organ <- "blood"
bl_16_up$organ <- "blood"
bl_18_up$organ <- "blood"
bl_18_down$organ <- "blood"
bl_20_down$organ <- "blood"
bl_20_up$organ <- "blood"
bl_term_up$organ <- "blood"
bl_term_down$organ <- "blood"


mu_4_up$wpi <- 4
mu_4_down$wpi <- 4
mu_8_up$wpi <- 8
mu_8_down$wpi <- 8
mu_12_up$wpi <- 12
mu_12_down$wpi <- 12
mu_14_up$wpi <- 14
mu_14_down$wpi <- 14
mu_18_up$wpi <- 18
mu_18_down$wpi <- 18
mu_20_down$wpi <- 20
mu_20_up$wpi <- 20
mu_term_up$wpi <- 22
mu_term_down$wpi <- 22
mu_4_up$organ <- "muscle"
mu_4_down$organ <- "muscle"
mu_8_up$organ <- "muscle"
mu_8_down$organ <- "muscle"
mu_12_up$organ <- "muscle"
mu_12_down$organ <- "muscle"
mu_14_up$organ <- "muscle"
mu_14_down$organ <- "muscle"
mu_18_up$organ <- "muscle"
mu_18_down$organ <- "muscle"
mu_20_down$organ <- "muscle"
mu_20_up$organ <- "muscle"
mu_term_up$organ <- "muscle"
mu_term_down$organ <- "muscle"



br_4_up$wpi <- 4
br_8_up$wpi <- 8
br_8_down$wpi <- 8
br_12_up$wpi <- 12
br_12_down$wpi <- 12
br_14_up$wpi <- 14
br_14_down$wpi <- 14
br_16_down$wpi <- 16
br_16_up$wpi <- 16
br_18_up$wpi <- 18
br_18_down$wpi <- 18
br_20_down$wpi <- 20
br_20_up$wpi <- 20
br_term_up$wpi <- 22
br_term_down$wpi <- 22
br_4_up$organ <- "brain"
br_8_up$organ <- "brain"
br_8_down$organ <- "brain"
br_12_up$organ <- "brain"
br_12_down$organ <- "brain"
br_14_up$organ <- "brain"
br_14_down$organ <- "brain"
br_16_down$organ <- "brain"
br_16_up$organ <- "brain"
br_18_up$organ <- "brain"
br_18_down$organ <- "brain"
br_20_down$organ <- "brain"
br_20_up$organ <- "brain"
br_term_up$organ <- "brain"
br_term_down$organ <- "brain"
br_4_up$Expressed_4wpi_main <- NULL
br_8_up$Expressed_8wpi_main <- NULL
br_8_down$Expressed_8wpi_main <- NULL
br_12_up$Expressed_12wpi_main <- NULL
br_12_down$Expressed_12wpi_main <- NULL
br_14_up$Expressed_14wpi_main <- NULL
br_14_down$Expressed_14wpi_main <- NULL
br_16_down$Expressed_16wpi_main <- NULL
br_16_up$Expressed_16wpi_main <- NULL
br_16_down$Expressed_18wpi_main <- NULL
br_16_up$Expressed_18wpi_main <- NULL
br_18_up$Expressed_18wpi_main <- NULL
br_18_down$Expressed_18wpi_main <- NULL
br_20_down$Expressed_20wpi_main <- NULL
br_20_up$Expressed_20wpi_main <- NULL
br_term_up$Expressed_term_main <- NULL
br_term_down$Expressed_term_main <- NULL
colnames(br_4_up) <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_8_up)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_8_down)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_12_up)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_12_down)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_14_up)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_14_down)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_16_down)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_16_up)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_18_up)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_18_down)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_20_down)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_20_up)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_term_up)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")
colnames(br_term_down)  <- c("gene_name","log2.Ratio","pValue","fdr","wpi","organ")

# merge all timepoints and organs
deg_down <- do.call("rbind",list(br_8_down,br_12_down,br_14_down,br_16_down,br_18_down,br_20_down,br_term_down,
                  sp_4_down,sp_16_down,sp_term_down,
                  mu_4_down,mu_8_down,mu_12_down,mu_14_down,mu_18_down,mu_20_down,mu_term_down,
                  bl_4_down,bl_14_down,bl_16_down,bl_18_down,bl_20_down,bl_term_down))
deg_up <- do.call("rbind",list(br_4_up,br_8_up,br_12_up,br_14_up,br_16_up,br_18_up,br_20_up,br_term_up,
                               sp_16_up,sp_16_up,sp_term_up,
                               mu_4_up,mu_8_up,mu_12_up,mu_14_up,mu_18_up,mu_20_up,mu_term_up,
                               bl_4_up,bl_16_up,bl_18_up,bl_20_up,bl_term_up))



# write tables and graph for monotonic hits
hif3a <- deg_up[deg_up$gene_name == "Hif3a",]
plin4 <- deg_up[deg_up$gene_name == "Plin4",]
c4b <- deg_up[deg_up$gene_name == "C4b",]
Hspa1b <- deg_down[deg_down$gene_name =="Hspa1b",]


hif3a$norm.ratio <- 2^hif3a$log2.Ratio
plin4$norm.ratio <- 2^plin4$log2.Ratio
c4b$norm.ratio <- 2^c4b$log2.Ratio
Hspa1b$norm.ratio <- 2^abs(Hspa1b$log2.Ratio)

write.table(hif3a,file="up_DEG.csv",append=FALSE,dec=".",sep=";")
write.table(plin4,file="up_DEG.csv",append=TRUE,dec=".",sep=";")
write.table(c4b,file="up_DEG_C4b.csv",append=TRUE,dec=".",sep=";")
write.table(Hspa1b,file="down_DEG_Hspa1b.csv",append=TRUE,dec=".",sep=";")

# make graphs for Plin4, Hspa1b, Col4a2, Hif3a
ggplot(deg_up[deg_up$gene_name == "Hif3a",],aes(x = wpi,y=log2.Ratio)) +
  xlim(0,22)+ylim(0.5,1.5) +
  geom_point(aes(color=organ),size=5)+
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"), 
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text=element_text(size=20),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))


ggplot(deg_up[deg_up$gene_name == "Plin4",],aes(x = wpi,y=log2.Ratio)) +
  xlim(0,22)+ylim(0.5,2.1) +
  geom_point(aes(color=organ),size=5)+
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"), 
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text=element_text(size=20),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))

ggplot(deg_down[deg_down$gene_name == "Hspa1b",],aes(x = wpi,y=log2.Ratio)) +
  xlim(0,22)+ylim(-2,0) +
  geom_point(aes(color=organ),size=5)+
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"), 
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text=element_text(size=20),
        legend.position = c(0.95,.8),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
        
ggplot(deg_down[deg_down$gene_name == "Col4a2",],aes(x = wpi,y=log2.Ratio)) +
  xlim(0,22)+ylim(-2,0) +
  geom_point(aes(color=organ),size=5)+
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"), 
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text=element_text(size=20),
        legend.position = c(0.95,.05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
