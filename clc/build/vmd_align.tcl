package require Orient
namespace import Orient::orient

# taken from here: http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/
set sel [atomselect top "all"]
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 1] {0 1 0}]
$sel move $A
set I [draw principalaxes $sel]

# additional rotation
set matrix [transaxis x 90]
$sel move $matrix

# save new pdb
$sel writepdb $argv
