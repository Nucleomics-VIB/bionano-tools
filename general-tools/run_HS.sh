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
# copy results to folder and create tar.gz archive; 2016/08/21; v1.1
#
# visit our Git: https://github.com/BITS-VIB

# check and adapt the following parameters for your system
# TIP: the paths may start with "/home/mic_common/" (with aliases in "/home/bionano")
# make sure you point to the latest code version !!
TOOLS="/home/bionano/tools"
SCRIPTS="/home/bionano/scripts"

#########################################
# please do not modify below this limit #
#########################################

version="1.0, 2016_10_22"

usage='# Usage: run_HS.sh
# script version '${version}'
## input files
# [required: -i <assembly-folder> (used for other paths below)]
# [required: -n <sequence fasta file>]
# [required: -b <BioNano CMAP file: exp_refineFinal1_contigs.cmap>]
# [required: -m <molecule BNX file to align molecules to genome maps and hybrid scaffolds>]
## conflict filtering: 1=no filter, 2=cut contig at conflict, 3=exclude conflicting contig 
# [-B <1|2|3 (filter for optical maps: default=2)>]
# [-N <1|2|3 (filter for sequences: default=2)>]
## required config settings with default values
# [-q <optArgument.xml (default to <assembly-folder>/optArguments_XXX.xml)>]
# [-e <errbin file (defaults to <assembly-folder>/output/contigs/auto_noise/autoNoise1.errbin)>]
# [-c <hybridScaffold_config.xml (default to <assembly-folder>/hybridScaffold_config.xml)>]
## other parameters with default values
# [-o <output folder (default to <assembly-folder>/hybridscaffold#>]
# [-p <path to Scripts (default to $SCRIPTS)>]
# [-s <hybridScafffold.pl file (default to $SCRIPTS/HybridScaffold/hybridScafffold.pl)>]
# [-r <RefAligner binary file (default to $TOOLS/RefAligner)>]
# [-h for this help]'

while getopts "i:n:b:m:B:N:q:e:o:c:p:s:r:h" opt; do
  case $opt in
    i) denovopath=${OPTARG} ;;
    n) fastaseq=${OPTARG} ;;
    b) refcmap=${OPTARG} ;;
	m) bnxfile=${OPTARG} ;;
    B) filtbnx=${OPTARG} ;;
    N) filtseq=${OPTARG} ;;
    q) optargpath=${OPTARG} ;;
    e) errbinfile=${OPTARG} ;;
    o) outpath=${OPTARG} ;;
    c) hybscafxml=${OPTARG} ;;
    p) scriptpath=${OPTARG} ;;
    s) hybridscapath=${OPTARG} ;;
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
	echo "! # ${1} file not found, provide with $2"
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
#errbin_path=${errbin_path:-"${denovo_path}/hybridscaffold/output/auto_noise/contigs/auto_noise/autoNoise1.errbin"}
errbin_path=${errbin_path:-$(find . -name "autoNoise1.errbin" -print | head -n 1 | sed -e 's/\.\///')}
testfileexist "${errbin_path}" "-e"

# check hybridScaffold_config.xml
#hybscaf_xml=${hybscafxml:-"${denovo_path}/hybridscaffold/hybridScaffold_config.xml"}
hybscaf_xml=${hybscafxml:-$(find . -name "hybridScaffold_config.xml" -print | head -n 1 | sed -e 's/\.\///')}
testfileexist "${hybscaf_xml}" "-c"

# check optArguments_XXX.xml
#optarg_path=${optargpath:-$(ls ${denovo_path}/hybridscaffold/optArguments*.xml)}
optarg_path=${optargpath:-$(find . -name "optArguments*.xml" -print | head -n 1 | sed -e 's/\.\///')}
testfileexist "${optarg_path}" "-q"

######################################
# optional settings exist when defined
######################################

filt_bnx=${filtbnx:-2}
filt_seq=${filtseq:-2}

##############################################
# create numbered hybridscaffold output folder
##############################################

out_path=${outpath:-"${denovopath}/hybridscaffold"}
if [[ -e "$out_path" ]] ; then
	i=2
	while [[ -e "$out_path-$i" ]] ; do
		let i++
	done
	name="$out_path-$i"
	out_path=${name}
fi
mkdir -f -p "$out_path"

# build command and quote weird chars
echo "# computing HS using the command:" | tee ${out_path}_log.txt

cmd="perl ${hybridscaff_path} \
	-n ${fasta_seq} \
	-b ${ref_cmap} \
	-m ${bnx_file} \
	-c ${hybscaf_xml} \
	-r ${refali_path} \
	-o ${out_path} \
	-f \
	-B ${filt_bnx} \
	-N ${filt_seq} \
	-x \
	-y \
	-p ${SCRIPTS} \
	-q ${optarg_path} \
	-e ${errbin_path} \
	2>&1 | tee -a ${out_path}_log.txt"
	
# store cmd in log
echo "# ${cmd}" 2>&1 | tee -a ${out_path}_log.txt

# execute cmd
${cmd}

if [ ! $? -eq 0 ]; then
    echo "! hybridscaffold command failed, please check your parameters"
	exit 1
fi

###############
# post process
###############

endts=$(date +%s)
dur=$(echo "${endts}-${startts}" | bc)
echo "HS Done in ${dur} sec" | tee -a ${out_path}_log.txt

echo "# now copying and archiving results"

# create result folder
base_folder=$(dirname ${out_path})
hs_base=$(basename ${out_path})
ref_base=$(basename ${ref_cmap%.cmap})
seq_base=$(basename ${fasta_seq%.f*})
result_folder=${base_folder}/hs_${hs_base}/${ref_base}_vs_${seq_base}_B${filt_bnx}_N${filt_seq}
mkdir -p ${result_folder}

# copy various input files
cp ${fasta_seq} ${result_folder}/
cp ${ref_cmap} ${result_folder}/
cp ${bnx_file} ${result_folder}/
cp ${hybscaf_xml} ${result_folder}/
cp ${optarg_path} ${result_folder}/
cp ${errbin_path} ${result_folder}/

# copy HS results
cp ${out_path}_log.txt ${result_folder}/
cp -r ${out_path}/* ${result_folder}/

# create archive from folder then delete original
tar -zcvf ${denovo_path}/${result_folder}.tgz ${result_folder} && rm -rf ${result_folder}
echo "# HS data was archived in ${result_folder}.tgz"

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
      