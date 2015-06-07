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

    # delete random existing ter cards
    sed -i '.bak' '/TER/d' $i
    # delete preface
    sed -i '.bak' '/TITLE/d' $i
    sed -i '.bak' '/REMARK/d' $i
    # delete model info
    sed -i '.bak' '/MODEL/d' $i
    sed -i '.bak' 's/ENDMDL/END/g' $i

    # ter cards to end of protein chains
    sed -i '.bak' 's/OT1/O  /g' $i
    sed -i '.bak' 's/OT2/OXT/g' $i
    sed -i '.bak' '/OXT/a\
    TER\
    '  $i

    # properly specify ions
    sed -i '.bak' 's/CLA/Cl-/g' $i
    sed -i '.bak' 's/POT/K+ /g' $i
    # ion ter cards
    sed -i '.bak' '/Cl- Cl-/a\
    TER\
    '  $i
    sed -i '.bak' '/K+  K+/a\
    TER\
    '  $i

    # water
    sed -i '.bak' 's/TIP3/WAT /g' $i
    sed -i '.bak' 's/OH2 WAT/O   WAT/g' $i
    sed -i '.bak' '/O   WAT/a\
    TER\
    '  $i

    # membrane
    # lipids
    sed -i '.bak' '/C118 OL/a\
    TER\
    '  $i
    sed -i '.bak' '/O12 PC/a\
    TER\
    '  $i
    sed -i '.bak' '/O12 PE/a\
    TER\
    '  $i
    sed -i '.bak' '/C116 PA/a\
    TER\
    '  $i

    # charmmgui artifact
    sed -i '.bak' 's/CD  ILE/CD1 ILE/g' $i

    # hist protonation
    sed -i '.bak' 's/HSD/HIE/g' $i
done
