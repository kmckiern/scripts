ffmpeg
	-i full_trj_macrocl_xy.mp4 -i full_trj_macrocl_z.mp4
	-filter_complex "
		nullsrc=size=1000x2250 [base];
		[0:v] setpts=PTS-STARTPTS, scale=1000x1125 [upperleft];
		[1:v] setpts=PTS-STARTPTS, scale=1000x1125 [upperright];
		[base][upperleft] overlay=shortest=1 [tmp1];
		[tmp1][upperright] overlay=shortest=1:x=1000
	"
	-c:v libx264 macro2.mp4
