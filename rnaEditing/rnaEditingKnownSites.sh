#!/bin/bash
# detect known editing sites using REDItoolKnown.py
SAMPLES='*.bam'
REF=/home/ubuntu/mm10/mm10.fa
KNOWN=/mnt/vol2/rnaediting/REDIportal/mm10_rediportal.txt.gz

for eachsample in $SAMPLES
do
        REDItoolKnown.py -i $eachsample -f $REF -l $KNOWN -o _REDItoolKnown -m 255 -t 14
done
