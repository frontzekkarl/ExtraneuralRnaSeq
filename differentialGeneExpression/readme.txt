Differentially expressed genes of extra-neural organs and brain from this and our previous publication https://doi.org/10.1371/journal.ppat.1008653, respectively, can be browsed freely via: XXXXXX.

Assessing differential gene expression was described previously here https://doi.org/10.1371/journal.ppat.1008653:

Quality control of reads was performed using FastQC. Low-quality ends were clipped (5’: 3 bases; 3’: 10 bases). Trimmed reads were aligned to the reference genome and transcriptome (FASTA and GTF files, respectively, downloaded from the UCSC mm10) with STAR version 2.3.0e_r291 with default settings.

Differentially expressed genes were identified based |log2FC| > 0 and FDR < 0.05 using the R package edgeR from Bioconductor (version 3.0). Only genes with at least 10 counts in at least 50% of the samples in one of the groups were considered in the analysis. Differentially expressed genes (DEGs) were defined as genes changing with |log2FC| > 0.5 and FDR < 0.05. An additional covariate was added to the model in edgeR to account for the batch effect associated with the different runs. 

R Session Info lists packages and versions used for DEG analysis.


R Scripts used for analysis are grouped as follows:
#############
#Main Cohort#
#############

main cohort: young mice inoculated with RML6 prions or non-infectious brain homogenate, organs harvested at 4,8,12,14,16,18,20 weeks post inoculation and terminal stage

blood -> bloodAnalysisMainCohort.R
muscle -> muscleAnalysisMainCohort.R
spleen -> spleenAnalysisMainCohort.R
