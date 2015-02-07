#!/bin/bash

for i in *jpg; do PREF=$(pref $i); convert -resize '800x800' -gravity center -background '#000000' -extent '800x800' $i ${PREF}_crop.jpg; done
