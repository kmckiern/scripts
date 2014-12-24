#!/bin/bash

# parse fasta sequences using a list of pdb names
# inpdb: list of input pdbs
# seqs: list of fasta sequences
# parsed: output file

INPDB=$1
SEQS=$2
PARSED=$3

while read LINE; do
    if [[ "$LINE" == ";"* ]]; then
        continue
    else
        grep -i -A 1 $LINE $SEQS >> $PARSED
    fi
done < $INPDB
