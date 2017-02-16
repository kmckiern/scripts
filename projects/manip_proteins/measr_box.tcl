set sel [atomselect top all]
#cell basis vectors
set m [measure minmax $sel]
echo $m "minmax results"
foreach {j1 j2} $m {}
foreach {x2 y2 z2} $j2 {}
foreach {x1 y1 z1} $j1 {}
set x_val [expr $x2 - $x1]
set y_val [expr $y2 - $y1]
set z_val [expr $z2 - $z1]
#cell origin
set origin [measure center $sel]
echo "xval" $x_val
echo "yval" $y_val
echo "zval" $z_val
echo "origin" $origin
quit
