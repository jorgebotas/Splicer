#!/bin/bash


coordinates=$1
refseq_genename=$2

echo $coordinates
echo $refseq_genename


join -1 2 -2 4 <(tail -n +2 $refseq_genename | sort -k 2,2) \
               <(sort -k 4,4 $coordinates) | head
