#!/bin/bash

# script name: run_SV.sh
# run BioNanoGenomics runSV.py with some help and testing.
# This bash script feeds parameters to the BNG python script runSV.py
# * it adds a number of additional tests to make sure all inputs are present
# * it facilitates the typing by resolving file paths automatically
# * it remains fully customisable like the original code

## Requirements:
# a unix computer installed with working bionano code and scripts (version > 2.5; 2.1)
# python present and working
# a folder with a full de-novo assembly file structure (from IrysSolve)
#
# Stephane Plaisance (VIB-NC+BITS) 2016/08/18; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# check and adapt the following parameters for your system
# TIP: the paths may start with "/home/mic_common/" (with aliases in "/home/bionano")
# make sure you point to the latest code version !!

# edit the following variables to match your system
#TOOLS="/home/bionano/tools"
#SCRIPTS="/home/bionano/scripts"

# try auto-detect (RefAligner is in PATH!)
TOOLS=$(dirname $(which RefAligner))
SCRIPTS=$(echo $TOOLS | sed -e 's/tools$/scripts/')

#########################################
# please do not modify below this limit #
#########################################

version="1.0, 2016_08_16"

usage='# Usage: run_SV.sh -r <path to the reference.cmap> -i <assembly-folder>
# script version '${version}'
# [optional: -t <path to Refaligner (default to "'$TOOLS'")>]
# [optional: -q <query_path containing assembly cmaps (default to <assembly>/output/contigs/exp_refineFinal1)>]
# [optional: -o <output folder (defaults in input map folder with same name + suffix "_sv")>]
# [optional: -p <path to the python script "runSV.py" (default to "'$SCRIPTS'")>]
# [optional: -a <optArgument.xml (default to "<assembly-folder>/optArguments_XXX.xml")>]
# [optional: -T <number of threads to use (default=4)>]
# [optional: -j <threads per job (default=4)>]
# [optional: -b <BED file with GAPs (unset by default)>]
# [optional: -e <.err file (defaults to <query_path>/alignref_final/EXP_REFINEFINAL1.err)>]
# [optional: -E <.errbin file (defaults to <query_path>/alignref_final/EXP_REFINEFINAL1.errbin)>]
# [optional: -C <CXML file for running on cluster (unset by default)>]
# [optional: -s <SV job configuration 0="single job" (default), 1="single job per contig" (not recommended), 2="grouped">
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

# execution time
startts=$(date +%s)

##########################################
# define default values for input settings
##########################################

ref_cmap=${refcmap}
denovo_path=${denovopath}
query_path=${querypath:-${denovo_path}/output/contigs/exp_refineFinal1}
err_path=${errpath:-${query_path}/alignref_final/EXP_REFINEFINAL1.err}
errbin_path=${errbinpath:-${query_path}/alignref_final/EXP_REFINEFINAL1.errbin}

# locate optArguments_XXX.xml in <denovo_path>
optarg_path=${optargpath:-$(ls ${denovo_path}/optArguments*.xml)}

################################
# optArguments_XXX.xml was found
################################

if [ -z "${optarg_path+x}" ]
then
   echo "# no optArguments*.xml found!"
   echo "${usage}"
   exit 1
fi

############################################
# server settings (adapt to your environment
############################################

refali_path=${refalipath:-$TOOLS}
if [[ ! -x "$refali_path/RefAligner" ]]
then
	echo "# RefAligner is not executable or found"
fi

script_path=${scriptpath:-$SCRIPTS}
if [ ! -f "${script_path}/runSV.py" ]
then
   	echo "# runSV.py file not found!";
   	echo "${usage}"
   	exit 1
fi

maxthr=${threadnum:-4}
thrperjob=${threadperjob:-4}

###########################################
# minimal arguments provided and path exist
###########################################

# ref_cmap file provided and exists
if [ -z "${ref_cmap+x}" ]
then
   echo "# no reference cmap provided!"
   echo "${usage}"
   exit 1
else
	if [ ! -f "${refcmap}" ]; then
    	echo "# ${refcmap} file not found!";
    	exit 1
    fi
fi

# denovo_path folder provided and exists
if [ -z "${denovo_path+x}" ]
then
	echo "# <assembly path> (containing the \'output\' folder)!"
	echo "${usage}"
	exit 1
else
	if [ ! -d "${denovo_path}" ]; then
   		echo "# ${denovo_path} folder not found!";
	    exit 1
	fi
fi

# query_path folder exists
if [ ! -d "${query_path}" ]; then
    echo "# ${query_path} folder not found!";
    exit 1
fi

# optarg_path file exists
if [ ! -f "${optarg_path}" ]; then
    echo "# ${optarg_path} file not found!";
    exit 1
fi

# err_path file exists
if [ ! -f "${err_path}" ]; then
    echo "# ${err_path} file not found!";
    exit 1
fi

# errbin_path file exists
if [ ! -f "${errbin_path}" ]; then
    echo "# ${errbin_path} file not found!";
    exit 1
fi

# ref_cmap provided and exists
if [ -z "${ref_cmap+x}" ]
then
   echo "# no reference cmap provided!"
   echo "${usage}"
   exit 1
fi

######################################
# optional settings exist when defined
######################################

# add bed if file provided and exists
if [ ! -z "${bedpath+x}" ]; then
	if [ ! -f "${bedpath}" ]; then
    	echo "# ${bedpath} file not found!";
    	exit 1
    else
    	bed_path="-b ${bedpath}"
	fi
else
	bed_path=''
fi

# add cxml if file provided and exists
if [ -z "$cxmlpath+x}" ]; then
	if [ ! -f "${cxmlpath}" ]; then
    	echo "# ${cxmlpath} file not found!";
    	exit 1
    else
    	cxml_path="-C ${cxmlpath}"
    fi
else
	cxml_path=''
fi

# add sv conf if file provided and exists
zero=0
if [ -z "$svconf+x}" ]; then
	if [[ $svconf -ge 0 && $svconf -le 2 ]] ; then
    	sv_conf="${svconf}"
	else
  		echo "# Invalid -s input: $svconf"
  		echo "${usage}"
  		exit 1
	fi
else
	sv_conf="${zero}"
fi

# create numbered output folder
bng_base=$(basename ${denovo_path})
ref_base=$(basename ${ref_cmap%.cmap})
out_path="${outpath:-${denovopath}/SV_${denovopath}_vs_${ref_base}}"
if [[ -e "$out_path" ]] ; then
	i=2
	while [[ -e "$out_path-$i" ]] ; do
		let i++
	done
	name=$out_path-$i
	out_path=${name}
fi

mkdir -p "${out_path}"

# from here down, redirect all outputs to log file
log_file="${out_path}/SV-analysis_log.txt"
touch ${log_file}

# build command and quote weird chars
echo "# computing SV from the provided data" | tee -a ${log_file}
echo | tee -a ${log_file}

cmd="python ${script_path}/runSV.py \
	-r ${ref_cmap} \
	-t ${refali_path} \
	-q ${query_path} \
	-o ${out_path} \
	-s ${sv_conf} \
	-p ${script_path} \
	-a ${optarg_path} \
	-T ${maxthr} \
	-j ${thrperjob} \
	-e ${err_path} \
	-E ${errbin_path} \
	${bed_path} \
	${cxml_path}"

echo "# ${cmd}" | tee -a ${log_file}

# execute cmd
{ ${cmd}; } 2>&1 | tee -a ${log_file}

if [ $? -ne 0 ] ; then
	echo "! SV-analysis command failed, please check your parameters" | \
		tee -a ${log_file}
	exit 1
fi

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)

echo | tee -a ${log_file}
echo "Done in ${dur} sec" | tee -a ${log_file}

###############
# post process
###############

echo "# now archiving results"
echo

# create archive from ${out_path} folder
arch_base=$(basename ${out_path})

# archive with tar and pigz if present
if hash pigz 2>/dev/null
then
	tar --use-compress-program="pigz -p8" -cvf ${denovopath}/${arch_base}.tgz ${out_path}
else
	tar -zcvf ${denovopath}/${arch_base}.tgz ${out_path}
fi

echo
echo "# MQR data was archived in ${denovopath}/${arch_base}.tgz"

exit 0

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
