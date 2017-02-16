#!/bin/bash

# for some reason amber requires ter cards in its pdbs
# however, vmd writes pdbs without ter cards
# this reinserts the ter cards
# it may not be perfect as flags were determined by eye

STDIN=( ${@} )

for i in "${STDIN[@]}"; do
    cp $i og.${i}
    # remove hydrogens
    /Applications/VMD\ 1.9.1.app/Contents/MacOS/startup.command -m $i -e ~/Dropbox/scripts/manip_proteins/rm_hydrogen.tcl -args $i

    # ter cards to end of protein chains
    sed -i '.bak' 's/OT1/O  /g' $i
    sed -i '.bak' 's/OT2/OXT/g' $i
    sed -i '.bak' '/OXT/a\
    TER\
    '  $i

    # charmmgui artifact
    sed -i '.bak' 's/CD  ILE/CD1 ILE/g' $i

    # hist protonation
    sed -i '.bak' 's/HSD/HIE/g' $i

    # lipids
    sed -i '.bak' '/C118 OL/a\
    TER\
    '  $i

    # ion ter cards
    sed -i '.bak' '/Cl- Cl-/a\
    TER\
    '  $i
    sed -i '.bak' '/K+  K+/a\
    TER\
    '  $i

    # water ter cards
    sed -i '.bak' '/O   WAT/a\
    TER\
    '  $i
done
