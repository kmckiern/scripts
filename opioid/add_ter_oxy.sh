#!/bin/bash

# sed -i '.bak' '/HH33 NME/a \
# TER\
# '  $1

sed -i '.bak' '/OXT ILE/a \
TER\
'  $1

sed -i '.bak' '/HO1 CHL/a \
TER\
'  $1

sed -i '.bak' '/H18T OL/a \
TER\
'  $1

# for NAL pdbs
sed -i '.bak' '/H23 NAL/a \
TER\
'  $1

# for OXY pdbs
sed -i '.bak' '/O1  OXY/a \
TER\
'  $1

sed -i '.bak' '/Cl- Cl-/a \
TER\
'  $1

sed -i '.bak' '/H2  WAT/a \
TER\
'  $1
