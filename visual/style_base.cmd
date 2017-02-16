windowsize 1000 1125

preset apply publication 1
light mode ambient

set shadows
set silhouette
set silhouetteWidth 3.8
set depthCue
set dcStart .4
set dcEnd .8
set flatTransparency

rep bs
setattr m stickScale .2
setattr m ballScale .12

sel
color byhetero sel

sel :Cl-
rep sphere sel
~show sel

setattr p drawMode 1
setattr p radius .06
setattr p color green
