#!/bin/bash

cat $1 | while read LINE; do
    LIG=$(echo $LINE | awk -F '/' '{print $1}')
    PREF=$(echo $LINE | awk -F '/' '{print $2}' | awk -F '.' '{print $1}')
    mkdir -p best/$PREF
    # i think individual dirs are need bc of temp files w non-unique names
    cp $LINE production/$PREF/$LIG.rst
done
