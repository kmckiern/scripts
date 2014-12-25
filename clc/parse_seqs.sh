#!/bin/bash

# parse fasta sequences for modeller using a list of pdb names
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
        echo ">P1;${LINE}" >> $PARSED
        echo "sequence:::::::::" >> $PARSED
        # specific to homodimers but whatevs
        grep -i -A 1 ${LINE}_A $SEQS | tail -n 1 >> $PARSED
        grep -i -A 1 ${LINE}_B $SEQS | tail -n 1 >> $PARSED
    fi
done < $INPDB
