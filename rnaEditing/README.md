
# RNA editing changes in prion-diseased mice

- for RNA editing analysis of individual transcripts, a list mm10_rediportal.txt.gz from REDIportal http://srv00.recas.ba.infn.it/atlas/ needs to be downloaded prior to analysis
- Fastq files from Kanata et al. can be found under GEO accession number GSE90977
- Analysis of pre-computed AEI and individual RNA transcript editing from REDItools can be found in sourceData/
- --- if provided AEI and REDItools output is used, skip first two steps of each instruction below

## How to re-analyze and replicate data from the manuscript

### A-to-I editing index

* Align only unique transcripts from fastq files: starAlignRnaEditing.sh
* Calculate A-to-I editing index (AEI): rnaEditingAei.sh
* Analysis of AEI for main cohorts of peripheral organs and Kanata et al. data: AEImainValidationKanataCohorts.ipynb 
* Analysis of AEI for brain from Sorce et al. https://doi.org/10.1371/journal.ppat.1008653: AeiPrionBrain.ipynb

### RNA editing of individual transcripts
* Align only unique transcripts from fastq files: starAlignRnaEditing.sh
* Generate list of RNA editing sites of individual RNA transcripts using REDItools: rnaEditingKnownSites.sh
* Filter editing sites according to Gabay 10.1038/s41467-022-28841-4: prionREDI.ipynb
* -------- after filtering, editing output needs to be exported to R for statistical analysis of differentially edited sites using REDIT-LLR
* -------- REDIT_LLR can be downloaded here https://github.com/gxiaolab/REDITs
* -------- script adapted to current manuscript, based on REDIT: AnalysisREDIT.R
* -------- REDIT output is further analysed in prionREDI.ipynb
