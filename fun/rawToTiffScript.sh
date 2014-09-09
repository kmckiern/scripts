#!/bin/sh

# Defaults:
# Album name.
NAME="Photos_$(date +"%Y_%m_%d")"
# Maximum dimension size of the resized image.
MAXD=1200
# Initial image format.
FORMAT0="CR2"
# Final image format.
FORMATF="tiff"
# P0: path to images on camera memory card.
P0="/Volumes/NONAME/DCIM/100CANON"
# .gif delay.
DELAY=100
GFLAG=0
# Help?
HFLAG=0

while getopts ":n:d:i:f:o:g:h" opt; 
do
	case $opt in
		n) 
			NAME=$OPTARG;;
		d) 
			MAXD=$OPTARG;; 
		i) 
			FORMAT0=$OPTARG;; 
		f) 
			FORMATF=$OPTARG;;
		o)
			P0=$OPTARG;;	
		g)
			GFLAG=1
			DELAY=$OPTARG;;
		h)
			HFLAG=1
			echo "
			-[flag]: [argument]. For all but h.
			-n: album name. 
			-d: maximum dimension for resized images. 
			-i: initial file format. 
			-f: final file format. 
			-o: original image directory.
			-g: gif delay (in ms).
		       	-h: flag descriptions.
			" >&2
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1;;
		:)
			echo "Missing argument for -$OPTARG" >&2
			exit 1;;
		esac
done

if [ $HFLAG -eq 0 ]
then
	# PBASE: common root among computer paths.
	PBASE="$HOME/Pictures/ABC/$NAME"
	# P1: path to the RAW directory on computer.
	P1="$PBASE/$FORMAT0"
	# P2: path to TIFF directory on computer.
	P2="$PBASE/$FORMATF"
	
	# Create file hierarchy by image type.
	mkdir -p $P1
	mkdir -p $P2
	
	# Move DCIM files to FORMAT0 directory.
	mv $P0/*.$FORMAT0 $P1/
	
	# Change to FORMATF directory, convert images from .FORMAT0 TO .FORMATF.
	cd $P1
	for i in *.$FORMAT0
	do
		IMAGENAME=$(echo $i | cut -d . -f 1)
		sips -s format $FORMATF -Z $MAXD $i --out $P2/$IMAGENAME.$FORMATF
	done
	
	gif() {
		# P3: path to GIF on computer.
		P3="$PBASE/gif"
		mkdir -p $P3
		cd $PBASE
		convert -delay $DELAY $P2/*.$FORMATF $P3/$NAME.gif
	}
	
	if [ $GFLAG -eq 1 ]
		then gif
	fi
fi
