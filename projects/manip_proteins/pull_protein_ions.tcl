# cut water and lipids
set sel [atomselect top "protein or ions"]
# write
$sel writepdb $argv
quit
