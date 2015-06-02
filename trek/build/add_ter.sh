#!/bin/bash

# for some reason amber requires ter cards in its pdbs
# however, vmd writes pdbs without ter cards
# this reinserts the ter cards
# it may not be perfect as flags were determined by eye

STDIN=( ${@} )

for i in "${STDIN[@]}"; do
    # delete random existing ter cards
    sed -i '.bak' '/TER/d' $i
    # delete preface
    sed -i '.bak' '/TITLE/d' $i
    sed -i '.bak' '/REMARK/d' $i
    # delete model info
    sed -i '.bak' '/MODEL/d' $i
    sed -i '.bak' 's/ENDMDL/END/g' $i

    # ter cards to end of protein chains
    sed -i '.bak' '/OT2 THR/a\
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
    sed -i '.bak' '/H2  WAT/a\
    TER\
    '  $i
done
