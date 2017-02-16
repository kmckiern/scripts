~sel

movie record supersample 2

wait 100

turn x -1 90
wait 90
freeze
wait 10

turn y 1 90
wait 90
freeze
wait 10

sel protein
surface probeRadius 2 sel
wait 50
transparency 100,s
wait 50

sel #0:84,85,86,87,88,126,127,128,129,130,380,381,382,383,384,476
color byhetero sel
~ribbon ~sel
wait 25
show sel
ribbackbone sel
wait 25

scale 1.005 100
wait 100

turn y 1 360
wait 360
freeze

wait 100

movie stop
movie encode format mp4 quality high output rotations.mp4

