#!/bin/bash

# unpack fah output zip dirs
while read TRJ; do tar xvfj $TRJ; done < ../lasts.dat
