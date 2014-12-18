#!/bin/bash

# hack a few pdb things

STDIN=( ${@} )

for i in "${STDIN[@]}"; do
    # end of protein 
    sed -i '.bak1' 's/ OT1 / O   /g' $i
    sed -i '.bak1' 's/ OT2 / OXT /g' $i

    # hist protonation
    sed -i '.bak1' 's/ HSD / HIS /g' $i

    # charmmgui artifact
    sed -i '.bak1' 's/CD  ILE/CD1 ILE/g' $i
done
