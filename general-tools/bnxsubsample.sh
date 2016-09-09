#!/bin/bash

# BioNanoGenomics BNX subsample and filter at CLI
#
# script name: bnx_subsample.sh
# Requirements:
# run on a unix computer installed with working bionano code

# Stephane Plaisance (VIB-NC+BITS) 2016/09/09; v1.0
#
# visit our Git: https://github.com/BITS-VIB

# check parameters for your system
TOOLS=$BNGTOOLS
version="1.0, 2016_09_09"

usage='# Usage: bnx_subsample.sh -i <molecules.bnx>
# script version '${version}'
# [optional: -o <output_prefix|bnx_subset>]
# [optional: -l <minlen|100>]
# [optional: -x <maxlen|2500>]
# [optional: -f <minlabels|5>]
# [optional: -g <maxlabels|200>]
# [optional: -a <maxai|0.6>]
# [optional: -s <minSNR|3.5>]
# [optional: -t <max-threads|24>]
# [optional: -m <max-ram|64>]
# [optional: -n <sample N molecules>]'

while getopts "i:o:l:x:f:g:a:s:t:m:n:h" opt; do
  case $opt in
    i) bnxdata=${OPTARG} ;;
    o) outprefix=${OPTARG} ;;
    l) minimumlen=${OPTARG} ;;
    x) maximumlen=${OPTARG} ;;
    f) minimumlabels=${OPTARG} ;;
    g) maximumlabels=${OPTARG} ;;
    a) maxaverageint=${OPTARG} ;;
    s) minmumsnr=${OPTARG} ;;
    t) maxthreads=${OPTARG} ;;
    m) maxmemory=${OPTARG} ;;
    n) molnumber=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
timestamp=$(date +%s)

# defaults
startts=$(date +%s)

# user-provided variables or defaults
minlen=${minimumlen:-100}
maxlen=${maximumlen:-2500}
minsites=${minimumlabels:-5}
maxsites=${maximumlabels:-200}
maxai=${maxaverageint:-0.6}
minsnr=${minmumsnr:-3.5}
maxthr=${maxthreads:-24}
maxmem=${maxmemory:-64}

# test if minimal arguments were provided
if [ -z "${bnxdata}" ]
then
   echo "# no bnx provided!"
   echo "${usage}"s
   exit 1
fi

if [ ! -f "${bnxdata}" ]; then
	echo "${bnxdata} file not found!"
	exit 1
fi

if [ -z "${molnumber}" ]
then
	lim=''
	outfile=${outprefix:-bnx_subset}
else
	lim="-randomize 1 -subset 1 ${molnumber}"
	outfile=${outprefix:-${molnumber}_bnx_subset}
fi

# build command and quote weird chars
echo "# creating subset from ${bnxdata}"
cmd="${TOOLS}/RefAligner \
	-i $(printf '%q' "${bnxdata}") \
	-o $(printf '%q' "${outfile}") \
	-merge \
	-bnx \
	-minlen ${minlen} \
	-maxlen ${maxlen} \
	-minsites ${minsites} \
	-maxsites ${maxsites} \
	-minSNR ${minsnr} \
	-MaxIntensity ${maxai} \
	-maxmem ${maxmem} \
	-maxthreads ${maxthr} \
	-stdout \
	-stderr \
	${lim}"

echo "# ${cmd}"
eval ${cmd}

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)
echo "Done in ${dur} sec"
