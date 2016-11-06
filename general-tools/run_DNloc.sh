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

# edit the following variables to match your system
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

# check RefAligner in <TOOLS>
refali_path="${TOOLS}/RefAligner"
testexecutable "${refali_path}"

# check Assembler in <TOOLS>
assembl_path="${TOOLS}/Assembler"
testexecutable "${assembl_path}"

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

# usage: pipelineCL.py [-h] [-T T] [-j MAXTHREADS] [-N N] [-G BED] [-i ITER]
#                      [-I IMG] [-b BNX] [-l LOCAL] [-t TOOLS] [-B BYPASS]
#                      [-e EXP] [-r REF] [-s GROUPSV] [-n] [-x] [-c CLEANUP]
#                      [-g SIZE] [-C CXML] [-w] [-a XML] [-L LAMBDAREF]
#                      [-p PERF] [-d] [-f NGSARG] [-F] [-u] [-U [GROUPCONTIGS]]
#                      [-v [VERSION]] [-S] [-V RUNSV] [-A] [-y] [-Y] [-m] [-H]
#
# Pipeline for de novo assembly - BioNano Genomics
#
# optional arguments:
#   -h, --help         show this help message and exit
#   -T T               Available threads per Node [default 1]
#   -j MAXTHREADS      Threads per job [default 1]
#   -N N               Number of split bnx files; number of pairwise jobs is
#                      N*(N-1)/2 (optional, default: -T)
#   -G BED             Bed file for gaps, used in structural variation (SV)
#                      detection to check for SV overlap with reference gaps
#   -i ITER            Number of extension and merge iterations (default=1, must
#                      be in range [0,10], use 0 to skip)
#   -I IMG             File with listed paths for image processing; no longer
#                      supported--do not use (use without .bnx)
#   -b BNX             Input molecule (.bnx) file, required
#   -l LOCAL           Location of output files root directory, required, will
#                      be created if does not exist; if does exist, will
#                      overwrite contents (may be error-prone)
#   -t TOOLS           Location of executable files (RefAligner and Assembler,
#                      required)
#   -B BYPASS          Skip steps, using previous result. <= 0:None,
#                      1:ImgDetect, 2:NoiseChar/Subsample, 3:Pairwise,
#                      4:Assembly, 5:RefineA, 6:RefineB, (7:RefineNGS -f only),
#                      7:merge0, 8+(i-1)*2:Ext(i), 9+(i-1)*2:Mrg(i),
#                      N+1:alignmol
#   -e EXP             Output file prefix (optional, default = exp)
#   -r REF             Reference file (must be .cmap), to compare resulting
#                      contigs (optional)
#   -s GROUPSV         SV jobs configuration: 0 = single job (required for
#                      correct haplotype calls), 1 = single job per contig (not
#                      recommended), 2 = grouped (default 0; optional)
#   -n                 Evaluate single molecule noise characterization
#   -x                 Exit after auto noise (noise characterization), do not
#                      preform de novo assembly
#   -c CLEANUP         Remove contig results (0 - keep all (default), 1 - remove
#                      intermediate files, 2 - store in sqlite, 3 - store in
#                      sqlite and remove)
#   -g SIZE            Organism genome size estimate in megabases, used for
#                      tuning assembly parameters [optional, if > 0, will modify
#                      parameters, if == 0, ignored, must be float]
#   -C CXML            Run on cluster, read XML file for submission arguments
#                      (optional--will not use cluster submission if absent)
#   -w                 Wipe clean previous contig results
#   -a XML             Read XML file for parameters (required)
#   -L LAMBDAREF       Lambda phage is spiked in, used for molecule scaling
#                      (only used with -I input)
#   -p PERF            Log performance in pipelineReport 0=None, 1=time, 2=perf,
#                      3=time&perf (default=1)
#   -d                 Retired option (contig subdirectories), always enabled.
#   -f NGSARG          Directory which contains NGS contigs as cmaps for use in
#                      refineNGS (no longer supported).
#   -F                 When using -f, disable min contigs check for all stages
#                      prior to refineNGS (no longer supported).
#   -u                 Do not perform final refinement (not recommended).
#   -U [GROUPCONTIGS]  Group contigs in refinement and extension stages [default
#                      ON, use 0 to disable]
#   -v [VERSION]       Print version; exit if argument > 1 supplied.
#   -S                 Log stdout/stderr of all RefAligner calls in files with
#                      same names as output files and suffix .stdout (default on
#                      --use this flag to turn off).
#   -V RUNSV           Detect structural variations. Default: only after final
#                      stage (normally refineFinal); if argument 2, also after
#                      refineB; if argument 0, disable.
#   -A                 Align molcules to final contigs (ON by default, use this
#                      to turn off).
#   -y                 Automatically determine noise parameters (requires
#                      reference; optional, default off)
#   -Y                 Disable scan scaling in auto noise (default on with
#                      reference and -y)
#   -m                 Disable molecule vs reference alignments (default on with
#                      reference)
#   -H                 Use HG19 (human genome) as reference, loaded from
#                      Analysis/SV/CopyNumberProfiles/DATA. Overrides -r
#                      argument.