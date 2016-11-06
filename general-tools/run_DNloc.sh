#!/bin/bash

# script name: run_DNloc.sh
# run BioNanoGenomics de-novo-assembly using only regular cpu's
# This bash script feeds parameters to $BNG_SCRIPTS/pipelineCL.py
# * it adds a number of additional tests to make sure all inputs are present
# * it facilitates the typing by resolving file paths automatically
# * it remains fully customisable like the original code

## Requirements:
# a unix computer installed with working bionano code and scripts (version > 2.5; 2.1)
# python present and working
#
# Stephane Plaisance (VIB-NC+BITS) 2016/08/18; v1.0
# add more control and parameter checking; 2016/10/25; v1.1
#
# visit our Git: https://github.com/BITS-VIB

# check and adapt the following parameters for your system
# TIP: the paths may start with "/home/mic_common/" (with aliases in "/home/bionano")
# make sure you point to the latest code version !!

TOOLS="/home/bionano/tools"
SCRIPTS="/home/bionano/scripts"
pipelineCL=$SCRIPTS/pipelineCL.py

#########################################
# please do not modify below this limit #
#########################################

version="1.0, 2016_11_06"

usage='# Usage: run_DNloc.sh
# script version '${version}'
## arguments
# [required: -b <molecule BNX file to assemble>]
# [required: -r <BioNano ref CMAP file for noise computation and stats>]
# [required: -x <optArgument.xml>]
# [optional: -o <assembly-base-folder (default current folder)>]
# [optional: -s <pipelineCL.py path (required if not in the default location)]
# [optional: -t <max-threads | 8>]
# [optional: -j <max-jobs (max-thread/2) | 4 >]
# [-h for this help]'

while getopts "b:r:x:o:s:t:j:h" opt; do
  case $opt in
    b) bnxpath=${OPTARG} ;;
    r) refcmappath=${OPTARG} ;;
    x) optargpath=${OPTARG} ;;
    o) denovopath=${OPTARG} ;;
    s) pipelineCLpath=${OPTARG} ;;
    t) maxthreads=${OPTARG} ;;
    j) maxjobs=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# execution time
startts=$(date +%s)

#########################
###### functions ########
#########################
function testvariabledef ()
{
if [ -z "${1}" ]
then
   echo "! # argument ${2} needs a value!"
   echo "${usage}"
   exit 1
else
	return 0
fi
}

function testfolderexist ()
{
if [ ! -d "${1}" ]
then
	echo "! # ${2} folder not found!"
	echo "${usage}"
	exit 1
else
	return 0
fi
}

function testfileexist ()
{
if [ ! -e "${1}" ]
then
	echo "! # ${1} file not found, provide with ${2}"
	echo "${usage}"
	exit 1
else
	return 0
fi
}

function testexecutable ()
{
if [[ ! -x "$1" ]]
then
	echo "! # ${1} is not executable or absent"
	echo "${usage}"
	exit 1
else
	return 0
fi
}

#############################
# test executables and inputs
#############################

# locate pipelineCL in <SCRIPTS>
pipelineCL=${pipelineCLpath:-"${SCRIPTS}/pipelineCL.py"}
testfileexist "${pipelineCL}" "-s"

# system limits (set for laptop, adapt up to tot-thread-number -2)
max_thr=${maxthreads:-4}

# max half the previous argument (dual thread jobs from '-N')
max_jobs=${maxjobs:-2}

###########################################
# minimal arguments provided and path exist
###########################################

# check denovo_path
denovo_path=${denovopath}
testfolderexist "${denovo_path}" "-o"

# check molecules.bnx
bnx_file=${bnxpath}
testvariabledef "${bnx_file}" "-b"
testfileexist "${bnx_file}" "-b"

# check ref-cmap
ref_cmap=${refcmappath}
testvariabledef "${ref_cmap}" "-r"
testfileexist "${ref_cmap}" "-r"

# check optArguments_haplotype.xml
opt_args=${optargpath}
testvariabledef "${opt_args}" "-x"
testfileexist "${opt_args}" "-x"

##############################################
# create numbered hybridscaffold output folder
##############################################

out_path=${outpath:-"${denovopath}/denovo_assembly"}

if [[ -e "$out_path" ]]
then
	i=2
	while [[ -e "$out_path-$i" ]]
	do
		let i++
	done
	name="$out_path-$i"
	out_path=${name}
fi

mkdir -p "${out_path}/output"

# copy data to ${out_path} folder
cp ${cmap_ref} ${out_path}/
cp ${bnx_file} ${out_path}/
cp ${opt_args} ${out_path}/

# from here down, redirect all outputs to log file
log_file="${out_path}/denovo-assembly_log.txt"
touch ${log_file}

echo "# $(date)" | tee -a ${log_file}
echo "# using pipelineCL.py version: $(python ${pipelineCL} -v | \
	grep "Pipeline Version:")" | tee -a ${log_file}
echo | tee -a ${log_file}

echo "# computing denovo assembly using the command:" | tee -a ${log_file}
echo | tee -a ${log_file}

cmd="python ${pipelineCL} \
 	-U \
 	-d \
 	-T ${max_thr} \
 	-j ${max_job} \
 	-N 2 \
 	-i 5 \
 	-a $BNG_SETTINGS/${opt_args} \
 	-w \
 	-t $BNG_TOOLS/ \
 	-l $outfolder/output \
 	-b $outfolder/${bnx_file} \
 	-r $outfolder/${cmap_ref}"

# print cmd to log
echo "# ${cmd}" | tee -a ${log_file}

echo | tee -a ${log_file}
echo "## using assembly settings from ${opt_args}" | tee -a ${log_file}
echo "## using molecules from ${bnx_file}" | tee -a ${log_file}
echo "## using reference cmap ${cmap_ref}" | tee -a ${log_file}
echo "## using ${max_thr} thread for ${max_job} jobs" | tee -a ${log_file}
echo | tee -a ${log_file}

# execute cmd
{ ${cmd}; } 2>&1 | tee -a ${log_file}

if [ $? -ne 0 ] ; then
	echo "! denovo assembly command failed, please check your parameters" | \
		tee -a ${log_file}
	exit 1
fi

###############
# post process
###############

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)

echo | tee -a ${log_file}
echo "denovo assembly Done in ${dur} sec" | tee -a ${log_file}

echo | tee -a ${log_file}
echo "# now copying raw data and archiving results" | tee -a ${log_file}

# create archive from ${out_path} folder
ref_base=$(basename ${ref_cmap%.cmap})
seq_base=$(basename ${fasta_seq%.f*})
arch_file=${ref_base}_vs_${seq_base}_B${filt_bnx}_N${filt_seq}.tgz

# compress using pigz if present
if hash pigz 2>/dev/null
then
	tar_option='--use-compress-program="pigz -p8"'
else
	tar_option=''
fi

# add selected files to archive
outf=$outfolder/output

find ${outf} -maxdepth 1 -type f -print0 | \
	xargs -0 { tar ${tar_option} \
	--exclude='*.tar.gz' \
	--exclude='*_of_*.bnx' \
	--exclude='*.map' \
	--exclude='*_refined_*' \
	--exclude='*_group*' \
	--exclude='all_sorted.bnx' \
	-cf ${denovo_path}/${arch_file} ;}

tar ${tar_option} \
	--exclude='*.tar.gz' \
	--exclude='*.map' \
	--exclude='*_refined_*' \
	--exclude='*_group*' \
	--append \
	--file=${denovo_path}/${arch_file} \
	${outf}/ref \
	${outf}/contigs/alignmolvref/merge \
	${outf}/contigs/exp_refineFinal1 \
	${outf}/contigs/exp_refineFinal1_sv

# test if 'copynumber' folder is present
if [ -d "${outf}/contigs/alignmolvref/copynumber" ]
	tar ${tar_option} \
		--exclude='*.tar.gz' \
		--exclude='*.map' \
		--exclude='*_refined_*' \
		--exclude='*_group*' \
		--append \
		--file=${denovo_path}/${arch_file} \
		${outf}/contigs/alignmolvref/copynumber \
else
	echo "# no copy number data available" | tee -a ${log_file}
fi

tar ${tar_option} \
	--append \
	--file=${denovo_path}/${arch_file} \
	${outf}/contigs/exp_refineFinal1/alignmol/merge/EXP_REFINEFINAL1_merge.map

echo | tee -a ${log_file}
echo "# denovo-assembly data was archived in ${arch_file}" | tee -a ${log_file}

exit 0
