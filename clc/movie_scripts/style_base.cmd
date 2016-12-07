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

sel :Cl-
color byhetero sel
rep sphere sel
show sel
