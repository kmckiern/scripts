#!/bin/bash

sed -i '.bak' '/H2  WAT X/a \
TER\
'  $1
