#!/bin/bash

# for some reason amber requires ter cards in its pdbs
# however, vmd writes pdbs without ter cards
# this reinserts the ter cards
# it may not be perfect as flags were determined by eye

STDIN=( ${@} )

for i in "${STDIN[@]}"; do
    bash ~/Dropbox/scripts/trek/build/ter_trek.sh $i

    # protonate gating GLUs
    sed -i '.bak' 's/GLU X 128/GLH X 128/g' $i
    sed -i '.bak' 's/GLU X 614/GLH X 614/g' $i
done
