#!/bin/bash
# 
# Usage: AlignScanImage --a1064 anc-1064-file --a1548 anc-1548-file -S start-line -M aligned-mask-file [-g gap-open-threshold]
# 
# a bash shell script that takes input arguments and calls a matlab function 
# to align the input scanning image with gene sequence alignment algorithm
# and output a alignment mask to a tif file.
#

# a function to show the usage of this shell script
function Usage()
{
    cat << EOF

$0 --a1064 anc-1064-file --a1548 anc-1548-file -S start-line -M aligned-mask-file [-g gap-open-threshold]

a bash shell script that takes input arguments and calls a matlab function to align the input scanning image with gene sequence alignment algorithm and output a alignment mask to a tif file.

options:

    -a --a1064    data anc file of 1064 nm
    -A --a1548    data anc file of 1548 nm
    -S           the index of the start line to be aligned
    -M           the mask file
    -g --gapopen the gap penalty, if not given, default value is -1024^2*2
    -h --help    show help info

Zhan Li, zhanli86@bu.edu

EOF
    exit 1
}

# set a default value for the optional argument -g
GAPOPEN=-1024^2*2

# read the options
#TEMP=`getopt -o a:A:s:l:S:M:g: --l a1064:,a1548:,ns:,nl:,start:,mask:,gapopen:: -- "$@"`
TEMP=`getopt -o a:A:S:M:g:h --l a1064:,a1548:,start:,mask:,gapopen::,help -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables
while true ; do
    case "$1" in
	-a|--a1064)
	    case "$2" in 
		"") shift 2 ;;
		*) ANC1064=$2 ; shift 2 ;;
	    esac ;;
	-A|--a1548)
	    case "$2" in 
		"") shift 2 ;;
		*) ANC1548=$2 ; shift 2 ;;
	    esac ;;
	# -s|--ns)
	#     case "$2" in 
	# 	"") shift 2 ;;
	# 	*) NS=$2 ; shift 2 ;;
	#     esac ;;
	# -l|--nl)
	#     case "$2" in 
	# 	"") shift 2 ;;
	# 	*) NL=$2 ; shift 2 ;;
	#     esac ;;
	-S|--start)
	    case "$2" in 
		"") shift 2 ;;
		*) START=$2 ; shift 2 ;;
	    esac ;;
	-M|--mask)
	    case "$2" in 
		"") shift 2 ;;
		*) MASK=$2 ; shift 2 ;;
	    esac ;;
	-g|--gapopen)
	    case "$2" in 
		"") GAPOPEN=-1024^2*2 ; shift 2 ;;
		*) GAPOPEN=$2 ; shift 2 ;;
	    esac ;;
	-h|--help)
	    Usage
	    exit 0
	    ;;
	--) shift ; break ;;
	*) echo "Internal error!" ; exit 1 ;;
    esac
done

# check the input parameters before starting processing
if [ -z "$ANC1064" ] ||  [ -z "$ANC1548" ] ||  [ -z "START" ] ||  [ -z "$MASK" ]; then
    Usage
    exit 1
fi

# display the input parameters
echo anc1064: $ANC1064
echo anc1548: $ANC1548
#echo ns: $NS
#echo nl: $NL
echo gapopen: $GAPOPEN
echo startline: $START
echo maskfile: $MASK

BASEDIR=$(dirname $0)
# run matlab program
ML="/usr/local/bin/matlab -nodisplay -nojvm -singleCompThread -r "
$ML "datacubefile1064='$ANC1064'; datacubefile1548='$ANC1548'; gapopen=$GAPOPEN; startline=$START; AlignedMaskFile='$MASK'; run $BASEDIR/AlignScanImageOnGeo"