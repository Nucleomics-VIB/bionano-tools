#!/bin/bash

# run BioNanoGenomics MQR at CLI
#
# Requirements:
# run on a unix computer installed with working bionano code
# a folder with at least one data.bnx file
# a reference.cmap file
#
# this script can be used for a batch of BNX files like:
# for b in *.bnx; do
#	run_MQR.sh -i $b -r myreference.cmap;
# done
#
# Stephane Plaisance (VIB-NC+BITS) 2015/03/27; v1.1
# added quoting paths to avoid issues with spaces
# added more filtering options, 2015/09/20; v2.0
# updated v2.4 parameters and added non-unique alignments, 2016-06-23; v2.1
#
# visit our Git: https://github.com/BITS-VIB

########################################################################
# default MQR parameters recommended for human samples in IrysView v2.4.
# -nosplit 2 -BestRef 1 -biaswt 0 -Mfast 0 -FP 1.5 -FN 0.15
# -sf 0.2 -sd 0.0 -A 5 -outlier 1e-3 -outlierMax 40 -
# endoutlier 1e-4 -S -1000 -sr 0.03 -se 0.2 -MaxSF 0.25 -MaxSE
# 0.5 -resbias 4 64 -maxmem 64 -M 3 3 - minlen 150 -T 1e-11
# -maxthreads 32 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 3 -hash
# -hashdelta 10 -hashoffset 1 -hashmaxmem 64 -insertThreads 4
# -maptype 0 -PVres 2 -PVendoutlier -AlignRes 2.0 -rres 0.9 -
# resEstimate -ScanScaling 2 -RepeatMask 5 0.01 -RepeatRec 0.7
# 0.6 1.4 -maxEnd 50 –usecolor 1 -stdout – stderr –randomize
# –subset 1 5000

# check parameters for your system
TOOLS=$BNG_TOOLS
version="2.1, 2016_06_23"

usage='# Usage: runMQR.sh -i <molecules.bnx> -r <reference.cmap>
# script version '${version}'
# [optional: -l <minlen|150>]
# [optional: -x <maxlen|2000>]
# [optional: -a <maxai|0.6>]
# [optional: -s <minSNR|3.5>]
# [optional: -p <pval|1e-9>]
# [optional: -u <BestRef (best-only=1; more=0)|1>]
# [optional: -b <if -u 0, #bestalignments|1>]
# [optional: -t <max-threads|24> -m <max-ram|64>]
# [optional: -n <sample N molecules>]'

while getopts "i:r:o:l:x:f:g:a:s:p:u:b:t:m:n:h" opt; do
  case $opt in
    i) bnxdata=${OPTARG} ;;
    r) refcmap=${OPTARG} ;;
    o) outfolder=${OPTARG} ;;
    l) minimumlen=${OPTARG} ;;
    x) maximumlen=${OPTARG} ;;
    f) minimumlabels=${OPTARG} ;;
    g) maximumlabels=${OPTARG} ;;
    a) maxaverageint=${OPTARG} ;;
    s) minmumsnr=${OPTARG} ;;
    p) pvalue=${OPTARG} ;;
    u) besthit=${OPTARG} ;;
    b) maxalign=${OPTARG} ;;
    t) maxthreads=${OPTARG} ;;
    m) maxmemory=${OPTARG} ;;
    n) molnumber=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
startts=$(date +%s)

# user-provided variables or defaults
minlen=${minimumlen:-150}
maxlen=${maximumlen:-2000}
minsites=${minimumlabels:-5}
maxsites=${maximumlabels:-200}
maxai=${maxaverageint:-0.6}
minsnr=${minmumsnr:-3.5}
pval=${pvalue:-"1e-9"}
best=${besthit:-1}

# add max alignment when BestRef (-u) is set to 0
if [ $best == 1 ]; then
	maxali=''
else
	maxali="-bestalignments ${opt_b:-1}"
fi

maxthr=${maxthreads:-24}
maxmem=${maxmemory:-64}

# test if minimal arguments were provided
if [ -z "${bnxdata}" ]
then
   echo "# no bnx provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f ${bnxdata} ]; then
	echo "${bnxdata} file not found!"
	exit 1
fi

if [ -z "${refcmap}" ]
then
	echo "# no reference cmap provided"
	echo "${usage}"
	exit 1
fi

if [ ! -f ${refcmap} ]; then
    echo "${refcmap} file not found!";
    exit 1
fi

if [ -z "${molnumber}" ]
then
	lim=''
else
	lim="-randomize 1 -subset 1 ${molnumber}"
fi

# create result folder
outpath=${outfolder:-'MQR-results'}

if [[ -e "$outpath" ]] ; then
	i=2
	while [[ -e $outpath-$i ]] ; do
		let i++
	done
	name=$outpath-$i
	outpath=${name}
fi

mkdir "$outpath"

# build command and quote weird chars
echo "# computing MQR from basename(${bnxdata})"
cmd="${TOOLS}/RefAligner -f \
	-ref ${refcmap} \
	-i $(printf '%q' "${bnxdata}") \
	-o $(printf '%q' "${outpath}/MoleculeQualityReport") \
	-nosplit 2 \
	-BestRef ${best} ${maxali} \
	-biaswt 0 \
	-Mfast 0 \
	-FP 1.5 \
	-FN 0.15 \
	-sf 0.2 \
	-sd 0.0 \
	-A 5 \
	-outlier 1e-3 \
	-outlierMax 40 \
	-endoutlier 1e-4 \
	-S -1000 \
	-sr 0.03 \
	-se 0.2 \
	-MaxSF 0.25 \
	-MaxSE 0.5 \
	-resbias 4 64 \
	-maxmem ${maxmem} \
	-M 3 3 \
	-minlen ${minlen} \
	-maxlen ${maxlen} \
	-minsites ${minsites} \
	-maxsites ${maxsites} \
	-minSNR ${minsnr} \
	-MaxIntensity ${maxai} \
	-T ${pval} \
	-maxthreads ${maxthr} 
	-hashgen 5 3 2.4 1.5 0.05 5.0 1 1 3 \
	-hash \
	-hashdelta 10 \
	-hashoffset 1 \
	-hashmaxmem ${maxmem} \
	-insertThreads 4 \
	-maptype 0 \
	-PVres 2 \
	-PVendoutlier \
	-AlignRes 2.0 \
	-rres 0.9 \
	-resEstimate \
	-ScanScaling 2 \
	-RepeatMask 5 0.01 \
	-RepeatRec 0.7 0.6 1.4 \
	-maxEnd 50 \
	-usecolor 1 \
	-stdout \
	-stderr \
	${lim}"

echo "# ${cmd}"
eval ${cmd}

endts=$(date +%s)
dur=$(echo ${endts}-${startts} | bc)
echo "Done in ${dur} sec"
