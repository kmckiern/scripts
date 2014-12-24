#!/bin/bash

# parse fasta sequences using a list of pdb names
# inpdb: list of input pdbs
# seqs: list of fasta sequences

INPDB=$1
SEQS=$2
PARSED=$3

while read LINE; do
    grep -A 1 $LINE $SEQS >> $PARSED
done < $INPDB
