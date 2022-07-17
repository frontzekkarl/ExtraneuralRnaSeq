args = commandArgs(trailingOnly = TRUE)
# args should be a character vector with Region , wpi and cohort elements

library(dplyr)
library(SGSeq)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DEXSeq)
library(org.Mm.eg.db)

# select the right samples 
all_samples_info <- readRDS("p3506_allSamples_info_main_cohort.rds")

#filter sample info based on arguments supplied in CL
selected_samples_info <- all_samples_info %>% filter(Region == args[1],
                                                     wpi == args[2])



# make a df with BAM file names  
bam_names <- data.frame(sample_name = selected_samples_info$SampleID,
                        file_bam = paste(selected_samples_info$SampleID, ".bam", sep = ""))

# extract info from BAM files
bam_info <- getBamInfo(sample_info = bam_names,
                       cores = 20)

# predict splice junctions and exons for each sample
tx_features <- predictTxFeatures(sample_info = bam_info,
                                 min_overhang = NULL,
                                 verbose = FALSE,
                                 cores = 20)

# merge predictions from individual samples to obtain a common set of transcript features
tx_features_merged <- mergeTxFeatures(tx_features)

# process terminal exons 
tx_features_merged_processed <- processTerminalExons(tx_features_merged)

# convert transcript features to splice graph features
sg_features <- convertToSGFeatures(tx_features_merged_processed)

# annotate splice graph features

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
seqlevelsStyle(txdb) <- "UCSC"
txf_ucsc <- convertToTxFeatures(txdb)

sg_features_anno <- annotate(sg_features, txf_ucsc)

#count features 

sg_feature_counts <- getSGFeatureCounts(sample_info = bam_info,
                                        features = sg_features_anno, 
                                        min_anchor = 1,
                                        counts_only = FALSE,
                                        verbose = FALSE,
                                        cores = 20)

#count variants
sg_variant_counts <- analyzeVariants(object = sg_feature_counts,
                                     cores = 20)

# sg_variant_counts stores 2 counts for each variant (one for the 5' and one for the 3' end of the variant). Get a single count for each variant
sgv <- rowRanges(sg_variant_counts)
sg_variant_counts_new <- getSGVariantCounts(sgv, sample_info = bam_info, cores = 20)

##  prepare input for DEXSeq -->  requires per variant counts, uique identifiers for each variant and a variable indicating how variants are grouped by events
sg_variant_cm <- counts(sg_variant_counts_new)
vid <- variantID(sg_variant_counts_new)
eid <- eventID(sg_variant_counts_new)

#keep only variants with at least five counts in at least three samples 

rows_to_keep <- function(counts){
  if(length(which(counts >= 5)) <3){
    x <- FALSE
  }else{
    x <- TRUE
  }
  x
}

good_rows_index <- apply(X = sg_variant_cm, MARGIN = 1, FUN = rows_to_keep)
sg_variant_cm_filtered <- sg_variant_cm[good_rows_index, ]
vid_filtered <- vid[good_rows_index]
eid_filtered <- eid[good_rows_index]


# find events that have only one variant left and remove them.

good_events_index <- !eid_filtered %in% as.integer(names(which(table(eid_filtered) == 1)))

eid_filtered2 <- eid_filtered[good_events_index]
vid_filtered2 <- vid_filtered[good_events_index]
sg_variant_cm_filtered2 <- sg_variant_cm_filtered[good_events_index,]


# create a df with sample info
selected_samples_info2 <- data.frame(row.names = factor(selected_samples_info$SampleID, levels= colnames(sg_variant_cm_filtered2)),
                                     condition = factor(selected_samples_info$Treatment))

selected_samples_info2 <- selected_samples_info2[order(row.names(selected_samples_info2)), , drop = FALSE ]


#construct a DEXSeqDataSet object:

dxd <- DEXSeqDataSet(countData =  sg_variant_cm_filtered2,
                     sampleData = selected_samples_info2,
                     design = ~ sample + exon + condition:exon,
                     featureID = as.character(vid_filtered2),
                     groupID = as.character(eid_filtered2))

#normalization

dxd = estimateSizeFactors(dxd)

# dispersion estimation
dxd = estimateDispersions(dxd)

# test for differential usage
dxd = testForDEU(dxd)
dxd = estimateExonFoldChanges(dxd, fitExpToVar = "condition")

# extract differential expression results 
dxr1 = DEXSeqResults(dxd)

#make a nice df with results 
dxr2 <- data.frame(dxr1[c(1:10)]) %>%
  dplyr::rename(eventID = groupID, variantID = featureID)


# map variants and events to genes and variant types 
sgvc_df <- mcols(sg_variant_counts_new)
variant_type_vec <- sapply(sgvc_df$variantType, FUN = paste, simplify = TRUE, collapse = "_")
geneName_vec <- sapply(sgvc_df$geneName, FUN = paste, simplify = TRUE, collapse = "_")
EV_gene_info_df <- data.frame(geneID = sgvc_df$geneID,
                              eventID = sgvc_df$eventID,
                              variantID = sgvc_df$variantID,
                              geneName = geneName_vec,
                              variantType = variant_type_vec)

dxr3 <- merge(dxr2, EV_gene_info_df, by = c("eventID", "variantID"), all.y = FALSE)
dxr3$geneSymbol <- mapIds(org.Mm.eg.db, keys = dxr3$geneName, keytype = "ENTREZID", column = "SYMBOL")


# find significant changes 

dxr3_signif <- dxr3 %>% filter(padj < 0.05)

saveRDS(dxr3_signif, "dxr3_signif.rds")



