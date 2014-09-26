#!/bin/bash

cat $1 | while read LINE; do
    LIG=$(echo $LINE | awk -F '/' '{print $1}')
    echo $LIG
    ls -ltrh $LINE
    tail -n 3 $LINE
done
