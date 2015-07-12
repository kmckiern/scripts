# cut water and lipids
set sel [atomselect top "protein"]
# write
$sel writepdb $argv

quit
