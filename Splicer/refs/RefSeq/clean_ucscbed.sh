#!/bin/bash

input=$1
refseq_genename=$2

if [ ! $refseq_genename ]; then
  echo "Refseq to genename correspondance must be provided"
  exit 1
fi

# Remove weird chromosomes
awk 'length($1) < 6' $input | 
  sed 's/_/./' | # Change transcript name (only include refseq id)
  sed 's/_/\t/g' | 
  sed 's/\./_/' | 
  cut -f 1-4,6,12 |
  sed 's/\./\t/' | # Remove transcript version
  cut -f 1-4,6,7 | # Only keep relevant information

  # Add gene name information to BED file
  join -1 2 -2 4 <(tail -n +2 $refseq_genename | sort -k 2,2) \
                 <(sort -k 4,4 -) | 

  awk 'BEGIN{ OFS="\t" }{ print $3, $4, $5, $1, $6, $7, $2 }'
