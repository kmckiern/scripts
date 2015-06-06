#!/bin/bash

# for some reason amber requires ter cards in its pdbs
# however, vmd writes pdbs without ter cards
# this reinserts the ter cards
# it may not be perfect as flags were determined by eye

STDIN=( ${@} )

for i in "${STDIN[@]}"; do
    # remove xtal lipids/ligands/non-k ions
    sed -i '.bak' '/PC1/d' $i
    sed -i '.bak' '/40D/d' $i
    sed -i '.bak' '/ CD /d' $i
done
