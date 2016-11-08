#!/bin/bash

# script name: run_HS.sh
# run BioNanoGenomics hybridScaffold.pl with some help and testing.
# This bash script feeds parameters to the BNG perl script hybridScaffold.pl
# * it adds a number of additional tests to make sure all inputs are present
# * it facilitates the typing by resolving file paths automatically
# * it remains fully customisable like the original code

## Requirements:
# a unix computer installed with working bionano code and scripts (version > 2.5; 2.1)
# python present and working
# a folder with a full de-novo assembly file structure (from IrysSolve)
#
# Stephane Plaisance (VIB-NC+BITS) 2016/08/18; v1.0
# add more control and parameter checking; 2016/10/25; v1.1
#
# visit our Git: https://github.com/BITS-VIB

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

version="1.1, 2016_10_25"

usage='# Usage: run_HS.sh
# script version '${version}'
## input files
# [required: -i <assembly-folder> (containing the output folder)]
# [required: -n <sequence fasta file>]
# [required: -b <BioNano CMAP file: EXP_REFINEFINAL1.cmap>]
# [required: -m <molecule BNX file to align molecules to genome maps and hybrid scaffolds>]
## conflict filtering: 1=no filter, 2=cut contig at conflict, 3=exclude conflicting contig 
# [-B <1|2|3 (filter for optical maps: default=2)>]
# [-N <1|2|3 (filter for sequences: default=2)>]
## required config settings with default values
# [-q <optArgument.xml (default to $SCRIPTS/optArguments_haplotype.xml)>]
# [-a use the hybscaf.xml _aggressive_ version (default OFF)]
# [-e <errbin file (defaults to <assembly-folder>/output/contigs/auto_noise/autoNoise1.errbin)>]
## other parameters with default values
# [-o <output folder (default to <assembly-folder>/hybridscaffold#>]
# [-p <path to Scripts (default to $SCRIPTS)>]
# [-s <hybridScafffold.pl file (default to $SCRIPTS/HybridScaffold/hybridScafffold.pl)>]
# [-r <RefAligner binary file (default to $TOOLS/RefAligner)>]
## by-default parameters or arguments not accessible using this script
# [-f and -x are set by default and not modifiable using this script]
# [-M cannot be set here (run secondary HS with manually edited conflicts.txt)]
# [-h for this help]'

while getopts "i:n:b:m:B:N:q:e:o:p:s:c:r:ah" opt; do
  case $opt in
    i) denovopath=${OPTARG} ;;
    n) fastaseq=${OPTARG} ;;
    b) refcmap=${OPTARG} ;;
    m) bnxfile=${OPTARG} ;;
    B) filtbnx=${OPTARG} ;;
    N) filtseq=${OPTARG} ;;
    a) aggressive=${OPTARG} ;;
    q) optargpath=${OPTARG} ;;
    e) errbinfile=${OPTARG} ;;
    o) outpath=${OPTARG} ;;
    p) scriptpath=${OPTARG} ;;
    s) hybridscapath=${OPTARG} ;;
    c) hybscafxml=${OPTARG} ;;
    r) refalipath=${OPTARG} ;;
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

# locate SCRIPTS
script_path=${scriptpath:-"${SCRIPTS}"}
testfolderexist "${script_path}" "scripts"

# locate HybridScaffold.pl in <SCRIPTS>
hybridscaff_path=${hybridscapath:-"${SCRIPTS}/HybridScaffold/hybridScaffold.pl"}
testfileexist "${hybridscaff_path}" "-s"

# locate RefAligner in <TOOLS>
refali_path=${refalipath:-"${TOOLS}/RefAligner"}
testexecutable "${refali_path}"

###########################################
# minimal arguments provided and path exist
###########################################

# check denovo_path
denovo_path=${denovopath}
testfolderexist "${denovo_path}" "-i"

# check sequence
fasta_seq=${fastaseq}
testvariabledef "${fasta_seq}" "-n"
testfileexist "${fasta_seq}" "-n"

# check cmap
ref_cmap=${refcmap}
testvariabledef "${ref_cmap}" "-b"
testfileexist "${ref_cmap}" "-b"

# check molecules.bnx
bnx_file=${bnxfile}
testvariabledef "${bnx_file}" "-m"
testfileexist "${bnx_file}" "-m"

# check autoNoise1.errbin
errbin_path=${errbin_path:-$(find ${denovopath} -name "autoNoise1.errbin" -print | \
	head -n 1 | sed -e 's/\.\///')}
testfileexist "${errbin_path}" "-e"

# check hybridScaffold_config.xml or hybridScaffold_config_aggressive.xml
if [ -z "${aggressive+x}" ]
then
	hybscaf_xml=${hybscafxml:-$(find ${script_path}/HybridScaffold -name "hybridScaffold_config.xml" -print | \
		head -n 1 | sed -e 's/\.\///')}
else
	hybscaf_xml=${hybscafxml:-$(find ${script_path}/HybridScaffold -name "hybridScaffold_config_aggressive.xml" -print | \
		head -n 1 | sed -e 's/\.\///')}
fi

# check if hybridScaffold_config file is present
testfileexist "${hybscaf_xml}" "-c"

# check optArguments_haplotype.xml
optarg_path=${optargpath:-$(find ${denovopath}/ -name "optArguments*.xml" -print | \
	head -n 1 | sed -e 's/\.\///')}
testfileexist "${optarg_path}" "-q"

######################################
# optional settings exist when defined
######################################

filt_bnx=${filtbnx:-2}
filt_seq=${filtseq:-2}

# create (0-based) array for reporting
declare -a filters=(null '1=no filter' '2=cut contig at conflict' '3=exclude conflicting contig')

##############################################
# create numbered hybridscaffold output folder
##############################################

out_path=${outpath:-"${denovopath}/hybridscaffold"}

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

# from here down, redirect all outputs to log file
log_file="${out_path}/HYBRID_SCAFFOLD_log.txt"
touch ${log_file}

echo "# $(date)" | tee -a ${log_file}
echo "# using hybridScaffold.pl version: $(perl ${hybridscaff_path} -v)" | tee -a ${log_file}
echo | tee -a ${log_file}

echo "# computing HS using the command:" | tee -a ${log_file}
echo | tee -a ${log_file}

# build hs command
cmd="perl ${hybridscaff_path} \
	-n ${fasta_seq} \
	-b ${ref_cmap} \
	-m ${bnx_file} \
	-c ${hybscaf_xml} \
	-r ${refali_path} \
	-o ${out_path}/output \
	-f \
	-B ${filt_bnx} \
	-N ${filt_seq} \
	-x \
	-p ${SCRIPTS} \
	-q ${optarg_path} \
	-e ${errbin_path} "

# print cmd to log
echo "# ${cmd}" | tee -a ${log_file}

echo | tee -a ${log_file}
echo "## filtering optical maps with "${filters["${filt_bnx}"]} | tee -a ${log_file}
echo "## filtering NGS maps with "${filters["${filt_seq}"]} | tee -a ${log_file}
echo "## using assembly settings from ${optarg_path} (use -q to change)" | tee -a ${log_file}
echo "## using scaffolding settings from ${hybscaf_xml} (use -c to change)" | tee -a ${log_file}
echo | tee -a ${log_file}

# execute cmd
{ ${cmd}; } 2>&1 | tee -a ${log_file}

if [ $? -ne 0 ] ; then
	echo "! hybridscaffold command failed, please check your parameters" | tee -a ${log_file}
	exit 1
fi

###############
# post process
###############

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)

echo | tee -a ${log_file}
echo "HS Done in ${dur} sec" | tee -a ${log_file}

echo | tee -a ${log_file}
echo "# now copying raw data and archiving results" | tee -a ${log_file}

# add raw data and archive
cp ${fasta_seq} ${out_path}/
cp ${ref_cmap} ${out_path}/
cp ${bnx_file} ${out_path}/
cp ${hybscaf_xml} ${out_path}/
cp ${optarg_path} ${out_path}/
cp ${errbin_path} ${out_path}/

# create archive from ${out_path} folder
ref_base=$(basename ${ref_cmap%.cmap})
seq_base=$(basename ${fasta_seq%.f*})
arch_file=HS-${ref_base}_vs_${seq_base}_B${filt_bnx}_N${filt_seq}.tgz

# archive with tar and pigz if present
if hash pigz 2>/dev/null
then
	tar --use-compress-program="pigz -p8" -cvf ${denovo_path}/${arch_file} ${out_path}
else
	tar -zcvf ${denovo_path}/${arch_file} ${out_path}
fi

echo | tee -a ${log_file}
echo "# HS data was archived in ${arch_file}" | tee -a ${log_file}

exit 0

# perl /home/mic_common/scripts_fresh/HybridScaffold/hybridScaffold.pl -h
#
# Usage: perl hybridScaffold.pl <-h> 
#	<-n ngs_file> <-b bng_cmap_file> <-c hybrid_config_xml> <-o output_folder> 
#	<-B conflict_filter_level> <-N conflict_filter_level> <-f> 
#	<-m molecules_bnx> <-p de_novo_pipeline> <-q de_novo_xml> <-v> <-x> <-y> <-e noise_param>
#       -h    : This help message         
#       -n    : Input NGS FASTA or CMAP file [required]
#       -b    : Input BioNano CMAP  [required]
#       -c    : Merge configuration file [required]
#       -o    : Output folder [required]
#       -r    : RefAligner program [required]
#       -B    : conflict filter level: 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contig [required if not using -M option]
#       -N    : conflict filter level: 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contig [required if not using -M option]
#       -f    : Force output and overwrite any existing files
#       -x    : Flag to generate molecules to hybrid scaffold alignment and molecules to genome map alignment [optional]
#       -y    : Flag to generate chimeric quality score for the Input BioNano CMAP [optional]
#       -m    : Input BioNano molecules BNX [optional; only required for either the -x or -y option]
#       -p    : Input de novo assembly pipeline directory [optional; only required for -x option]
#       -q    : Input de novo assembly pipeline optArguments XML script [optional; only required for -x option]
#       -e    : Input de novo assembly noise parameter .errbin or .err file [optional; recommended for -y option but not required]
#       -v    : Print pipeline version information
#       -M    : Input a conflict resolution file indicating which NGS and BioNano conflicting contigs to be cut [optional]
