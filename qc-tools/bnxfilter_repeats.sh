#!/bin/bash

# run BioNanoGenomics repeat filtering at CLI
#
# script name: bnxfilter_repeats.sh
# Requirements:
# run on a unix computer installed with working bionano code

# Stephane Plaisance (VIB-NC+BITS) 2016/03/04; v1.0
#
# visit our Git: https://github.com/BITS-VIB

# check parameters for your system
TOOLS=/home/bionano/tools

usage='# Usage: bnxfilter_repeats.sh
#  -i <input (bnx file)>
## optional parameters (|default value)
#  -t <stretch tolerance|0.1>
#  -m <min Repeat Units|5>
#  -c <choice (1/2/3)|1>
#  -l <keep log>'

while getopts "i:t:m:c:lh" opt; do
  case $opt in
    i) input=${OPTARG} ;;
    t) tolerance=${OPTARG} ;;
    m) minele=${OPTARG} ;;
    c) choice=${OPTARG} ;;
    l) log=TRUE ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
timestamp=$(date +%s)

# user-provided variables or defaults
tol=${tolerance:-0.1}
min=${minele:-5}
act=${choice:-1}
if [ ${log} ]; then
	keeplog="-stderr -stdout"
fi

# test if minimal arguments were provided
if [ -z "${input}" ]
then
   echo "# no bnx provided!"
   echo "${usage}"
   exit 1
fi

# test input files
if [ ! -f ${input} ]; then
    echo "${input} file not found!";
    exit 1
else
	# deduced variables
	foldername=$(dirname "${input}")
	name=$(basename "${input}" .bnx)
fi

# manage action
case ${act} in
	[1]) tag="nonrepeats"; ;;
	[2]) tag="repeats"; ;;
	[3]) tag="repeatmask"; ;;
	*) echo "invalid value for -c (1,2,3)"; exit 1 ;;
esac;

# build command
#  -simpleRepeatTolerance <X>
#		: fraction which adjacent intervals are allowed to differ to be considered
#		the same repeat element [Default: 0.1]
#	-simpleRepeatMinEle <X>
#		: minimum number of repeat elements required [Default: 5]
#	-simpleRepeatStandalone
#		: same as -simpleRepeat, but all alignments are skipped
#		(and their output files are not produced), and -merge is disregarded.
#	-simpleRepeatFilter <1,2,3>
#		: if non-zero, output filtered maps:
#		1 = output only maps which do not contain any repeats
#		2 = output only maps which contain any repeats,
#		3 = output all maps but with all repeats masked (ie, labels removed)
#	code also outputs .rmap describing repeats found,
#   	and repeatStats.txt file summarising repeat statistics.

echo "## filtering ${name} for ${choice}:${tag}"

cmd="${TOOLS}/RefAligner \
	-i ${input} \
	-f -simpleRepeatFilter ${act} \
	-o ${foldername}/${name} \
	-simpleRepeatTolerance ${tol} \
	-simpleRepeatMinEle ${min} \
	-simpleRepeatStandalone ${keeplog}"

echo -e "## command: ${cmd}"
eval ${cmd}

# did the command succeed?
if [ $? -eq 0 ]; then
	echo "## results were saved to ${foldername}/${name}_${tag}.bnx"
	echo "## additional info was written to ${foldername}/${name}.rmap and ${foldername}/${name}_repeatStats.txt"
	echo
	if [ -f "${foldername}/${name}_repeatStats.txt" ]; then
		echo
		echo "## summary"
		grep -v "^#" ${foldername}/${name}_repeatStats.txt | column -s ":" -t
	fi
else
    echo "# the command failed!"
fi
