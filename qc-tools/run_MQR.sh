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
#
# visit our Git: https://github.com/BITS-VIB

# check parameters for your system
TOOLS=/home/bionano/tools

usage='# Usage: runMQR.sh -i <molecules.bnx> -r <reference.cmap>
#		[optional: -l <minlen|150> -p <pval|1e-9>]
#		[optional: -t <max-threads|32> -m <max-ram|64>]
#		[optional: -s <sample N molecules>]'

while getopts "i:r:l:p:t:m:s:h" opt; do
  case $opt in
    i) bnx=${OPTARG} ;;
    r) ref=${OPTARG} ;;
    l) minlen=${OPTARG} ;;
    p) pval=${OPTARG} ;;
    t) maxthr=${OPTARG} ;;
    m) maxmem=${OPTARG} ;;
    s) molnum=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
timestamp=$(date +%s)

# user-provided variables or defaults
minlen=${opt_l:-150}
pval=${opt_p:-"1e-9"}
maxmem=${opt_m:-64}
maxthr=${opt_t:-32}

# test if minimal arguments were provided
if [ -z "${bnx}" ]
then
   echo "# no bnx provided!"
   echo "${usage}"
   exit 1
fi

if [ -z "${ref}" ]
then
   echo "# no ref provided"
   echo "${usage}"
   exit 1
fi

if [ -z "${molnum}" ]
then
	lim=''
else
	lim="-randomize 1 -subset 1 ${molnum}"
fi

# test input files
if [ ! -f ${bnx} ]; then
    echo "${bnx} file not found!";
    exit 1
else
	# deduced variables
	foldername=$(dirname "${bnx}")
<<<<<<< HEAD
	name=$(basename ${bnx} .bnx)
=======
	name=$(basename "${bnx}" .bnx)
>>>>>>> origin/master
fi

if [ ! -f ${ref} ]; then
    echo "${ref} file not found!";
    exit 1
fi

# build command and quote weird chars
echo "# computing MQR from ${name}.bnx"
cmd="${TOOLS}/RefAligner -f \
	-ref ${ref} \
	-i $(printf '%q' "${bnx}") \
<<<<<<< HEAD
    -o $(printf '%q' "${foldername}/MoleculeQualityReport") \
=======
	-o $(printf '%q' "${foldername}/MoleculeQualityReport") \
>>>>>>> origin/master
	-nosplit 2 -BestRef 1 -biaswt 0 -Mfast 0 -FP 1.5 -sf 0.2 -sd 0.0 -A 5 \
	-outlier 1e-4 -endoutlier 1e-3 -S -1000 -sr 0.04 -resbias 5 64 \
	-maxmem ${maxmem} \
	-M 3 \
	-minlen ${minlen} \
	-T ${pval} \
	-maxthreads ${maxthr} \
	-hashgen 5 3 2.4 1.5 0.05 5.0 1 1 2 -hash -hashdelta 10 \
	-hashmaxmem ${maxmem} \
	-insertThreads 8 \
	-stdout -stderr \
	${lim}"

echo "# ${cmd}"
eval ${cmd}
