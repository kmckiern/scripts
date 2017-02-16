#!/bin/sh

# $1 is the output dir.
# $2 is trj start.
# $3 is trj end.
# $4 is the number of boots.

# Produces A_L(T, t) and K_{Al}(T, t), as well as averages.
# I'd like to also add something for uncertainty.

for j in 3*; do echo -ne "Box-X\nBox-Y\n\n" | g_energy -f lipid-md.edr -b $2 -e $3 -o $1/al_$j.xvg; done

for i in 3*; do python get_kappa_b.py $1/al_$i.xvg 128 $i $1 $4 >> $1/result.dat; done
