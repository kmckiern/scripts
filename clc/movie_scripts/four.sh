ffmpeg
	-i full_trj_macrocl_xy.mp4 -i full_trj_macrocl_x.mp4 -i tic_trj.mp4 -i full_trj_macrocl_z.mp4
	-filter_complex "
		nullsrc=size=2000x2250 [base];
		[0:v] setpts=PTS-STARTPTS, scale=1000x1125 [upperleft];
		[1:v] setpts=PTS-STARTPTS, scale=1000x1125 [upperright];
		[2:v] setpts=PTS-STARTPTS, scale=1000x1125 [lowerleft];
		[3:v] setpts=PTS-STARTPTS, scale=1000x1125 [lowerright];
		[base][upperleft] overlay=shortest=1 [tmp1];
		[tmp1][upperright] overlay=shortest=1:x=1000 [tmp2];
		[tmp2][lowerleft] overlay=shortest=1:y=1125 [tmp3];
		[tmp3][lowerright] overlay=shortest=1:x=1000:y=1125
	"
	-c:v libx264 output.mkv
