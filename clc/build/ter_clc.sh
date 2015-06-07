#!/bin/bash

# for some reason amber requires ter cards in its pdbs
# however, vmd writes pdbs without ter cards
# this reinserts the ter cards
# it may not be perfect as flags were determined by eye

STDIN=( ${@} )

for i in "${STDIN[@]}"; do
    cp ${i} wH_${i}
    # remove hydrogens
    /Applications/VMD\ 1.9.1.app/Contents/MacOS/startup.command -m wH_${i} -e ~/Dropbox/scripts/trek/build/rm_hydrogen.tcl -args ${i}

    # clean up pdb
    # delete random existing ter cards
    sed -i '.bak' '/TER/d' $i
    # delete preface
    sed -i '.bak' '/TITLE/d' $i
    sed -i '.bak' '/REMARK/d' $i
    # delete model info
    sed -i '.bak' '/MODEL/d' $i
    sed -i '.bak' 's/ENDMDL/END/g' $i

    # end of protein 
    sed -i '.bak' '/OT2 PRO/a\
    TER\
    '  $i

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
    sed -i '.bak' '/O12 PGR/a\
    TER\
    '  $i

    # ions
    sed -i '.bak' '/Cl- Cl-/a\
    TER\
    '  $i
    sed -i '.bak' '/CAL CAL/a\
    TER\
    '  $i
    sed -i '.bak' '/C0  C0/a\
    TER\
    '  $i
    sed -i '.bak' '/K+  K+/a\
    TER\
    '  $i
    
    # protonation state of histidine
    sed -i '.bak' 's/HSD/HIE/g' $i
    # protonate gating GLUs
    sed -i '.bak' 's/GLU X 128/GLH X 128/g' $i
    sed -i '.bak' 's/GLU X 614/GLH X 614/g' $i

    # because amber is weird
    sed -i '.bak' 's/CD  ILE/CD1 ILE/g' $i

    # terminal oxygens
    sed -i '.bak' 's/OT1/O  /g' $i
    sed -i '.bak' 's/OT2/OXT/g' $i
    
    # water
    sed -i '.bak' '/O   WAT/a\
    TER\
    '  $i
done
