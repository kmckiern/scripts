#!/bin/bash

# for some reason amber requires ter cards in its pdbs
# however, vmd writes pdbs without ter cards
# this reinserts the ter cards
# it may not be perfect as flags were determined by eye

STDIN=( ${@} )

for i in "${STDIN[@]}"; do
    # remove hydrogens
    /Applications/VMD\ 1.9.1.app/Contents/MacOS/startup.command -m $i -e ~/Dropbox/scripts/trek/build/rm_hydrogen.tcl -args noH_${i}

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
    sed -i '.bak' '/H18T OL/a\
    TER\
    '  $i

    sed -i '.bak' '/O12 PC/a\
    TER\
    '  $i

    sed -i '.bak' '/O12 PE/a\
    TER\
    '  $i

    sed -i '.bak' '/H16T PA/a\
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
    
    # water
    sed -i '.bak' '/H2  WAT/a\
    TER\
    '  $i
done
