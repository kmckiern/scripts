#!/bin/bash

STDIN=( ${@} )

for i in "${STDIN[@]}"; do
    # sed -i '.bak' '/HH33 NME/a \
    # TER\
    # '  $i
    
    sed -i '.bak' '/OXT ILE/a \
    TER\
    '  $i
    
    sed -i '.bak' '/HO1 CHL/a \
    TER\
    '  $i
    
    sed -i '.bak' '/H18T OL/a \
    TER\
    '  $i
    
    # for NAL pdbs
    sed -i '.bak' '/H23 NAL/a \
    TER\
    '  $i
    
    # for OXY pdbs
    sed -i '.bak' '/O1  OXY/a \
    TER\
    '  $i
    
    sed -i '.bak' '/Cl- Cl-/a \
    TER\
    '  $i
    
    sed -i '.bak' '/H2  WAT/a \
    TER\
    '  $i
done
