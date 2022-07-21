
# RNA editing changes in prion-diseased mice

- for RNA editing analysis of individual transcripts, a list mm10_rediportal.txt.gz from REDIportal http://srv00.recas.ba.infn.it/atlas/ needs to be downloaded prior to analysis

### How to re-analyze and replicate data from the manuscript

* Align fastq files from the manuscript and Kanata et al., GEO accession number GSE90977: starAlignRnaEditing.sh
* Calculate A-to-I editing index (AEI): rnaEditingAei.sh
* Analysis of AEI for main cohorts including brain and Kanata et al. data: AEImainValidationKanataCohorts.ipynb and AeiPrionBrain.ipynb
* Generate list of RNA editing of individual RNA transcripts using REDItools: rnaEditingKnownSites.sh
