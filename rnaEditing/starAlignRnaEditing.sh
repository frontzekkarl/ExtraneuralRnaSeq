#!/bin/bash

# --readfilesCommand,silence out for unzipped files, unmute for gzipped
# add this to STAR accordingly:
# --readFilesCommand zcat \

# read all samples from folder
SAMPLES='*R1.fastq.gz'

for eachsample in $SAMPLES
do
	BASE="${eachsample::-11}"
 	R1=${BASE}R1.fastq.gz
  R2=${BASE}R2.fastq.gz
  echo "dis R1: $R1"
	echo "dis R2: $R2"

for eachsample in $SAMPLES
do

  STAR \
	--runThreadN 16 \
	--genomeDir ~/mm10/ \ # reference genome
	--genomeLoad LoadAndKeep \
	--limitBAMsortRAM 1200000000000 \
	--outFileNamePrefix $BASE \
	--outReadsUnmapped Fastx \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMstrandField intronMotif \
	--outSAMattributes All \
	--readFilesCommand zcat \
	--outFilterType BySJout \
	--outFilterMultimapNmax 1 \ # unique alignments only
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
  --readFilesIn $R1 $R2 
done
#STAR --genomeLoad Remove # free up memory after alignment
