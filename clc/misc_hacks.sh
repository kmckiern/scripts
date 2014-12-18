#!/bin/bash

# hack a few pdb things

STDIN=( ${@} )

for i in "${STDIN[@]}"; do
    # end of protein 
    sed -i '.bak' 's/ OT1 / O   /g' $i
    sed -i '.bak' 's/ OT2 / OXT /g' $i

    # hist protonation
    sed -i '.bak' 's/ HSD / HIS /g' $i

    # charmmgui artifact
    sed -i '.bak' 's/CD  ILE/CD1 ILE/g' $i
done
