#!/bin/bash

# remove non-protein non-ion molecules from pdb

STDIN=( ${@} )

for i in "${STDIN[@]}"; do
    # remove xtal lipids/ligands/non-k ions
    sed -i '.bak' '/PC1/d' $i
    sed -i '.bak' '/40D/d' $i
    sed -i '.bak' '/ CD /d' $i
    sed -i '.bak' '/B40/d' $i
    sed -i '.bak' '/A40/d' $i
    sed -i '.bak' '/A405/d' $i
    sed -i '.bak' '/A408/d' $i
    sed -i '.bak' '/TRD/d' $i
done
