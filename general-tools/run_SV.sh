#!/bin/bash

# run BioNanoGenomics SV calling at CLI
#
# Requirements:
# run on a unix computer installed with working bionano code and scripts
# a project where a full de-novo assembly dataset is present (IrysSolve)
#
# Stephane Plaisance (VIB-NC+BITS) 2016/08/18; v1.0
#
# visit our Git: https://github.com/BITS-VIB

########################################################################
# usage: runSV.py [-h] [-t REFALIGNER] [-r REFERENCEMAP] [-q QUERYDIR]
#                 [-o OUTPUTDIR] [-p PIPELINEDIR] [-a OPTARGUMENTS]
#                 [-T NUMTHREADS] [-j MAXTHREADS] [-b BEDFILE] [-e ERRFILE]
#                 [-E ERRBINFILE] [-C CXML] [-s GROUPSV]
# 
# Standalone script for running the SV Module of the Pipeline, ie, aligning
# genome maps (ie, bioNano contigs as .cmap) to sequence contigs or a reference
# (.cmap) for the purpose of structural variation (SV) detection.
# 
# optional arguments:
#   -h, --help       show this help message and exit
#   -t REFALIGNER    Path to RefAligner or dir containing it (required)
#   -r REFERENCEMAP  Path to reference maps (.cmap), 1 file only (required)
#   -q QUERYDIR      Path to dir containing query maps (.cmaps) (required)
#   -o OUTPUTDIR     output dir (optional, defaults to input map dir with suffix
#                    "_sv")
#   -p PIPELINEDIR   Pipeline dir (optional, defaults to script dir, or current
#                    directory)
#   -a OPTARGUMENTS  Path to optArguments.xml (optional, default
#                    optArguments_human.xml in Pipeline dir if found, otherwise
#                    required)
#   -T NUMTHREADS    Total number of threads (cores) to use (optional, default
#                    4)
#   -j MAXTHREADS    Threads per Job, -maxthreads (non-cluster only;optional,
#                    default 4)
#   -b BEDFILE       .bed file with gaps in reference for flagging SVs which
#                    overlap N-base gaps (optional)
#   -e ERRFILE       .err file to use for noise parameters--will supersede noise
#                    parameters in the optArgument supplied (but that file must
#                    still be supplied for non-noise parameters)
#   -E ERRBINFILE    .errbin file to use for noise parameters--will supersede
#                    noise parameters in the optArgument supplied (but that file
#                    must still be supplied for non-noise parameters)
#   -C CXML          Run on cluster, read XML file for submission arguments
#                    (optional--will not use cluster submission if absent)
#   -s GROUPSV       SV jobs configuration: 0 = single job (required for correct
#                    haplotype calls), 1 = single job per contig (not
#                    recommended), 2 = grouped (default 0; optional)


# check and adapt the following parameters for your system
TOOLS="/home/bionano/tools"
SCRIPTS="/home/bionano/scripts"

#########################################
# please do not modify below this limit #
#########################################

version="1.0, 2016_08_16"

usage='# Usage: runMQR.sh -r <path to the reference.cmap> -i <path to the de-novo assembly>
# script version '${version}'
# [optional: -t <path to Refaligner (opt:default to "/home/bionano/tools")>]
# [optional: -q <path to folder containing the query cmaps >]
# [optional: -o <output folder (defaults to input map dir with suffix "_sv")>]
# [optional: -p <path to the python script "runSV.py" (default "/home/bionano/scripts")>]
# [optional: -a <optArgument.xml (default /home/bionano/scripts/optArguments_human.xml)>]
# [optional: -T <number of threads to use (default 4)>]
# [optional: -j <threads per job (default 4)>]
# [optional: -b <BED file with GAPs (default none)>]
# [optional: -e <.err file>]
# [optional: -E <.errbin file>]
# [optional: -C <CXML file for running on cluster (default none)>]
# [optional: -s <SV job configuration 0=single job, 1="single job per contig (not recommended), 2 grouped (default 0)>
# [-h for this help]'

while getopts "r:i:t:q:o:p:a:T:j:b:e:E:C:s:h" opt; do
  case $opt in
    r) refcmap=${OPTARG} ;;
    i) denovopath=${OPTARG} ;;
    t) refalipath=${OPTARG} ;;
    q) querypath=${OPTARG} ;;
    o) outpath=${OPTARG} ;;
    p) scriptpath=${OPTARG} ;;
    a) optargpath=${OPTARG} ;;
    T) threadnum=${OPTARG} ;;
    j) threadperjob=${OPTARG} ;;
    b) bedpath=${OPTARG} ;;
    e) errpath=${OPTARG} ;;
    E) errbinpath=${OPTARG} ;;
    C) cxmlpath=${OPTARG} ;;
    s) svconf=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
startts=$(date +%s)

# input settings
ref_cmap=${refcmap}
denovo_path=${denovopath}
query_path=${querypath:-${denovo_path}/output/contigs/exp_refineFinal1}
err_path=${errpath:-${query_path}/alignref_final/EXP_REFINEFINAL1.err}
errbin_path=${errbinpath:-${query_path}/alignref_final/EXP_REFINEFINAL1.errbin}

optarg_path=${optargpath:-$SCRIPTS/optArguments_haplotype.xml}

# server settings
refali_path=${refalipath:-$BNGTOOLS}
script_path=${scriptpath:-$SCRIPTS}
maxthr=${threadnum:-4}
thrperjob=${threadperjob:-4}

# add bed if file provided
if [ ! -z "${bedpath}" ]; then
	bed_path="-b ${bedpath}"
else
	bed_path=''
fi

# add cxml if file provided
if [ -z "$cxmlpath}" ]; then
	cxml_path="-C ${cxmlpath}"
else
	cxml_path=''
fi

# add sv conf if file provided
if [ -z "$svconf}" ]; then
	sv_conf="-s ${svconf}"
else
	sv_conf=''
fi

# test if minimal arguments were provided
if [ -z "${ref_cmap}" ]
then
   echo "# no reference cmap provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${refcmap}" ]; then
    echo "${refcmap} file not found!";
    exit 1
fi

if [ -z "${denovo_path}" ]
then
   echo "# no denovo assembly path provided (containing the \'output\' folder)!"
   echo "${usage}"
   exit 1
fi

if [ ! -d "${query_path}" ]; then
    echo "${query_path} folder not found!";
    exit 1
fi

if [ ! -f "${err_path}" ]; then
    echo "${err_path} file not found!";
    exit 1
fi

if [ ! -f "${optarg_path}" ]; then
    echo "${optarg_path} file not found!";
    exit 1
fi

if [ ! -f "${err_path}" ]; then
    echo "${err_path} file not found!";
    exit 1
fi

if [ ! -f "${errbin_path}" ]; then
    echo "${errbin_path} file not found!";
    exit 1
fi

# output folder
out_path="${outpath:-${denovopath}/${denovopath}_sv}"

if [[ -e "$out_path" ]] ; then
	i=2
	while [[ -e "$out_path-$i" ]] ; do
		let i++
	done
	name=$out_path-$i
	out_path=${name}
fi

mkdir -p "$out_path"

# build command and quote weird chars
echo "# computing SV from the provided data"
cmd="python ${SCRIPTS}/run_SV.py \
	-r ${ref_cmap} \
	-t ${refali_path} \
	-q ${query_path} \
	-o ${out_path} \
	-p ${script_path} \
	-a ${optarg_path} \
	-T ${maxthr} \
	-j ${thrperjob} \
	${bed_path} \
	-e ${err_path} \
	-E ${errbin_path} \
	${cxml_path} \
	${sv_conf}"

echo "# ${cmd}"
eval ${cmd}

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)
echo "Done in ${dur} sec"
