set sel [atomselect top "not hydrogen"] 
$sel writepdb $argv
quit
