#!/bin/bash
# script to calculate A-I editing index according to Eisenberg et al, Nature Methods
# INBAM: folder of bams
INBAM=$1

RNAEditingIndex -d $INBAM/bam -l $INBAM/bam/AEI/log -o $INBAM/bam/AEI/output -os $INBAM/bam/AEI/summary -f _HpAligned.sortedByCoord.out.bam --genome mm10
