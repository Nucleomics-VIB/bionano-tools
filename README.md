[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
![bionano-tools](pictures/Irys_icon.png) - BioNano-Tools
==========

*All tools presented below have only been tested by me and may contain bugs, please let me know if you find some. Each tool relies on dependencies normally listed at the top of the code (cpan for perl and cran for R will help you add them)*

Please refer to the accompanying **[wiki](https://github.com/Nucleomics-VIB/bionano-tools/wiki)** for examples and workflows.

### Table of Contents
**[QC-tools](#qc-tools)**

- **[bnxstats.pl](#bnxstatspl)** - **[bnx2quantiles.pl](#bnx2quantilespl)** - **[bnx2quantiles.pl](#bnx2quantilespl)** - **[bnxfilter.pl](#bnxfilterpl)** - **[bnxfilter2.pl](#bnxfilter2pl)** - **[bnxsubsample.sh](#bnxsubsamplesh)** - **[bnxfilter_repeats.pl](#bnxfilter_repeatspl)** - **[bnxreheader.pl](#bnxreheaderpl)**  - **[run_mqr.sh](#run_mqrsh)**  

**[General-tools](#general-tools)**

- **[xmapisec.pl](#xmapisecpl)** - **[mapisec.pl](#mapisecpl)** - **[bnx0convert.pl](#bnx0convertpl)** - **[mqr2bnx.pl](#mqr2bnxpl)** - **[bnxclnheader.pl](#bnxclnheaderpl)** - **[bnxreheader.pl](#bnxreheaderpl-1)** - **[bedrename.pl](#bedrenamepl)** - **[cmap2bed.pl](#cmap2bedpl)** - **[xmap2bed.pl](#xmap2bedpl)** - **[xmap2bed12.pl](#xmap2bed12pl)** - **[smap2bed.pl](#smap2bedpl)** - **[bnxsplitter.pl](#bnxsplitterpl)** - **[labeldensity.pl](#labeldensitypl)** - **[cmap2renum.pl](#cmap2renumpl)** - **[findnregions.pl](#findnregionspl)** - **[fastafiltlength.pl](#fastafiltlengthpl)** - **[fastasortlength.pl](#fastasortlengthpl)** - **[fastarename.pl](#fastarenamepl)** - **[fastaPad20k.pl](#fastapad20kpl)** - **[FastaExtractList.pl](#fastaextractlistpl)** - **[CmapExtractList.pl](#cmapextractlistpl)** - **[run_SV.sh](#run_svsh)** - **[run_HS.sh](#run_hssh)** - **[run_DNloc.sh](#run_dnlocsh)**

**[SysAdmin-tools](#sysadmin-tools)**  

- **[logphicards.sh](#logphicardssh)** -

**[Shiny-apps](#shiny-apps)**  

- **[BNX_viewer.shinyapp](#bnx_viewershinyapp)** -


## QC-tools
*[[back-to-top](#top)]*  

**[bnxstats.pl](#bnxstatspl)** - **[bnx2quantiles.pl](#bnx2quantilespl)** - **[bnx2quantiles.pl](#bnx2quantilespl)** - **[bnxfilter.pl](#bnxfilterpl)** - **[bnxfilter2.pl](#bnxfilter2pl)** - **[bnxsubsample.sh](#bnxsubsamplesh)** - **[bnxfilter_repeats.pl](#bnxfilter_repeatspl)** - **[bnxreheader.pl](#bnxreheaderpl)**  - **[run_mqr.sh](#run_mqrsh)**  

Tools to process and QC BioNanoGenomics data.

*Please do dry run or argument '-h' to discover the available options in each command*

### **bnxstats.pl**
*[[QC-tools](#qc-tools)]*  

The perl script **[bnxstats_v1.pl](qc-tools/bnxstats_v1.pl)** computes basic stats from a BNX file with size filtering just like what IrysView does under windows (but without the W). This original script was renamed 'bnxstats_v1' and left for reference. A new version of the script **[bnxstats.pl](qc-tools/bnxstats.pl)** allows filtering by size but also by avgIntensity and SNR and returns counts for each filtering subset, it also supports gzipped BNX files which is a nice thing given the huge size these text files tend to develop.
```bash
## Usage: bnxstats.pl <-i bnx-file>
# Additional optional parameters are:
# <-l minsize in kb (100)>
# <-x maxsize in kb (5000)>
# <-s minimal SNR (3.50)>
# <-p percentile (99)>
# <-m max-AvgIntensity (or percentile if undef)>
# <-h to display this help>
```

### **bnx2quantiles.pl**
*[[QC-tools](#qc-tools)]*  

The perl script **[bnx2quantiles.pl](qc-tools/bnx2quantiles.pl)** analyze BNX data and return value distributions for the most important measurements. Useful to define cutoffs to be used with **bnxfilter.pl** and **bnxfilter2.pl**. 

```bash
You must provide a BNX file with -i
## Usage: bnx2quantiles.pl <-i bnx-file>
# script version:1.2
# Additional optional parameters are:
# <-p additional low percentile (1)>
# <-P additional high percentile (99)>
# <-h to display this help>
```
example run with BNG demo EColi data
```
# BioNanoGenomics data quantiles (v1.2), 06/13/2016
# input file: RawMolecules.bnx
# tot-molecules: 3'112
# inputFile size (Gb): 0.277
# additional percentiles:
#   low_percentile: 1
#   high_percentile: 99
#--------------------------------------------------

# distributions for All labelled molecules (N=3'099, 99.6%)
                   .        min        25%     median        75%        max       mean      stdev        N50         1%        99%   skewness   kurtosis
         length (kb)          1         33         70        121        626         89         76        130          3        358          0          0
   molecAvgIntensity       0.01       0.07       0.09       0.11       0.67       0.09       0.04       0.10       0.02       0.24       3.31      28.05
         molecAvgSNR       0.70       6.96      11.11      14.83      55.93      11.37       5.78      13.67       2.02      26.34       1.01       3.52
        labelDensity       1.34      29.18      33.30      37.70     109.80      33.71      10.10      34.82       8.83      67.72       1.14       6.59
    labelAvgDistance       0.00    2590.52    2930.23    3289.90   11062.50    2950.95     841.76    3043.53       0.00    5682.45       0.80       9.52
   labelAvgIntensity       0.00       0.02       0.03       0.05       0.51       0.04       0.03       0.05       0.00       0.15       3.04      25.55
         labelAvgSNR       0.32       4.65       7.27      10.31      42.16       7.78       4.82       9.74       0.84      22.41       1.33       4.63
#--------------------------------------------------
```

### **bnxfilter.pl**
*[[QC-tools](#qc-tools)]*  

The perl script **[bnxfilter.pl](qc-tools/bnxfilter.pl)** filters BNX data based on min- and max-length, max-averageIntensity, min-SNR to generate nicer data for assembly. The avgIntensitty value can be read from a run with **bnxstats.pl**, a default absolute value of '0.6' will otherwise be applied. The most recent version of the code supports gzipped data and exports as gzip as option.
```bash
## Usage: bnxfilter.pl <-i bnx-file>
# Additional optional parameters are:
# <-l min-length in kb (100)>
# <-x max-length in kb (5000)>
# <-m max-AvgIntensity (0.6)>
# <-s min-SNR (3.5)>
# <-n min-nicks (6)>
# <-z zip results (default OFF)>
```

### **bnxfilter2.pl**
*[[QC-tools](#qc-tools)]*  

The perl script **[bnxfilter2.pl](qc-tools/bnxfilter2.pl)** adds to the first version and allows filtering on **label average-intensity** and **label average-snr** (REM: both new parameters are optional and undef by default, please note the difference between molecule [upper-case] and label [lower-case] arguments).
```bash
You must provide a BNX file with -i
## Usage: bnxfilter2.pl <-i bnx-file>
# script version:2.0
# Additional optional parameters are:
# <-l min-length in kb (100)>
# <-x max-length in kb (2500)>
# <-S min-molSNR (3.5)>
# <-M max-molAvgIntensity (0.6)>
# <-n min-nicks (6)>
# <-s min-labSNR (undef)>
# <-m max-labAvgIntensity (undef)>
# <-z zip results (default OFF)>
# <-h to display this help>
```

### **bnxsubsample.sh**
*[[QC-tools](#qc-tools)]*  

The perl script **[bnxsubsample.sh](qc-tools/bnxsubsample.sh)** also filters molecules but makes use of the Bionano RefAligner executable instead of perl and is much faster; It also creates a subset (random) of molecules that can be usefull when less data is needed.

```bash
bnxsubsample.sh -i <molecules.bnx>
# script version 1.0, 2016_09_09
# [optional: -o <output_prefix|bnx_subset>]
# [optional: -l <minlen|100>]
# [optional: -x <maxlen|2500>]
# [optional: -f <minlabels|5>]
# [optional: -g <maxlabels|200>]
# [optional: -a <maxai|0.6>]
# [optional: -s <minSNR|3.5>]
# [optional: -t <max-threads|24>]
# [optional: -m <max-ram|64>]
# [optional: -n <sample N molecules>]
```

### **bnxfilter_repeats.pl**
*[[QC-tools](#qc-tools)]*  

The bash script **[bnxfilter_repeats.sh](qc-tools/bnxfilter_repeats.sh)** filters BNX data to 'remove', 'restrict to' or 'mask' simple repeats. It reflects the Windows version found in Irysview (that generates data with a wrongly formatted header) and works only on your linux server as it makes direct use of RefAligner. The code is simplistic and you could as well type the command in your terminal.
```bash
# Usage: bnxfilter_repeats.sh
# script version 2.1, 2016_06_23
#  -i <input (bnx file)>
## optional parameters (|default value)
#  -t <stretch tolerance|0.1>
#  -m <min Repeat Units|5>
#  -c <choice (1/2/3)|1>
#  -l <keep log>
#
# -c 1 = output only maps which do not contain any repeats
# -c 2 = output only maps which contain any repeats,
# -c 3 = output all maps but with all repeats masked (ie, labels removed)
```

### **bnxreheader.pl**
*[[QC-tools](#qc-tools)]*  

The perl script **[bnxreheader.pl](qc-tools/bnxreheader.pl)** replaces unsupported characters in the header of a BNX file by '_' to avoid issue in IrysView (eg MQR returning empty sample name when other chatacters are present). This script should become obsolete when BNG correct their code or validate user input.
```bash
## Usage: bnxfilter.pl <-i bnx-file>
# Additional optional parameters are:
# <-z zip results (default OFF)>
# <-h to display this help>
```

### **run_MQR.sh**
*[[QC-tools](#qc-tools)]*  

The bash script **[run_MQR.sh](qc-tools/run_MQR.sh)**

Perform molecule quality report (MQR) at CLI instead of running this under IrysView. One may prefer to perform the MQR directly on his/her Nix server. The main advantage is that one can launch this code in a bash loop and perform all MQR from a list of BNX files (single runs) without supervision.

Parameters added to the MQR comands match those found in the documentation for Irys version 2.4. Some of these parameters WILL change in future version, please check with your current Irys version.

```
# default MQR parameters recommended for human samples in IrysView v2.4.
-nosplit 2 -BestRef 1 -biaswt 0 -Mfast 0 -FP 1.5 -FN 0.15
-sf 0.2 -sd 0.0 -A 5 -outlier 1e-3 -outlierMax 40 -
endoutlier 1e-4 -S -1000 -sr 0.03 -se 0.2 -MaxSF 0.25 -MaxSE
0.5 -resbias 4 64 -maxmem 64 -M 3 3 - minlen 150 -T 1e-11
-maxthreads 32 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 3 -hash
-hashdelta 10 -hashoffset 1 -hashmaxmem 64 -insertThreads 4
-maptype 0 -PVres 2 -PVendoutlier -AlignRes 2.0 -rres 0.9 -
resEstimate -ScanScaling 2 -RepeatMask 5 0.01 -RepeatRec 0.7
0.6 1.4 -maxEnd 50 –usecolor 1 -stdout – stderr –randomize
–subset 1 5000
```

The following loop will process all BNX files in the current folder and create a MQR (in its own folder) from each BNX and a common reference (more parameters are available).
```bash
for b in *.bnx; do
  run_MQR.sh -i $b -r myreference.cmap;
done
```

Type the script name followed by -h will list all available parameters
```bash
# Usage: runMQR.sh -i <molecules.bnx> -r <reference.cmap>
# script version 2.1, 2016_06_23
# [optional: -o <outfolder (default to current folder)>]
# [optional: -l <minlen|150>]
# [optional: -x <maxlen|2000>]
# [optional: -a <maxai|0.6>]
# [optional: -s <minSNR|3.5>]
# [optional: -p <pval|1e-9>]
# [optional: -u <BestRef (best-only=1; more=0)|1>]
# [optional: -b <if -u 0, #bestalignments|1>]
# [optional: -t <max-threads|24>]
# [optional: -m <max-ram|64>]
# [optional: -n <sample N molecules>]
```

## General-tools
*[[back-to-top](#top)]*  

**[xmapisec.pl](#xmapisecpl)** - **[mapisec.pl](#mapisecpl)** - **[bnx0convert.pl](#bnx0convertpl)** - **[mqr2bnx.pl](#mqr2bnxpl)** - **[bnxclnheader.pl](#bnxclnheaderpl)** - **[bnxreheader.pl](#bnxreheaderpl-1)** - **[bedrename.pl](#bedrenamepl)** - **[cmap2bed.pl](#cmap2bedpl)** - **[xmap2bed.pl](#xmap2bedpl)** - **[xmap2bed12.pl](#xmap2bed12pl)** - **[smap2bed.pl](#smap2bedpl)** - **[bnxsplitter.pl](#bnxsplitterpl)** - **[labeldensity.pl](#labeldensitypl)** - **[cmap2renum.pl](#cmap2renumpl)** - **[findnregions.pl](#findnregionspl)** - **[fastafiltlength.pl](#fastafiltlengthpl)** - **[fastasortlength.pl](#fastasortlengthpl)** - **[fastarename.pl](#fastarenamepl)** - **[fastaPad20k.pl](#fastapad20kpl)** - **[FastaExtractList.pl](#fastaextractlistpl)** - **[FastaExtractList.pl](#fastaextractlistpl)** - **[CmapExtractList.pl](#cmapextractlistpl)** - **[run_SV.sh](#run_svsh)** - **[run_HS.sh](#run_hssh)** - **[run_DNloc.sh](#run_dnlocsh)**

### **xmapisec.pl**
*[[General-tools](#general-tools)]*  

The perl script **[xmapisec.pl](general-tools/xmapisec.pl)** takes information from two MQR runs (from the 'MoleculeQualityReport.xmap' file) to split the corresponding BNX file into molecules that align to either, both, or none of the reference cmaps used for eather MQR. This allows creating subset of a BNX file that may be more specific for one or another reference assembly (which could for example represent parental genomes for a diploid) and make the resulting BNX data accessible for other applications like denovo assembly. An optional parameter allows filtering alignments by their 'Alignment Score' to create more stringent datasets. The MQR runs may be performed using a lower than normal '-T' value in order to keep molecules that align with some degree of divergence (parental genomes are not necessariy identical to the haploid compound of a hybrid genome). *We would welcome your feedback after using this tool and reports of success would be a great reward for the work put into this script.*
```
Aim: Identify molecules specific to two ref-cmaps, ubiquitous, or not-aligning

# Usage: xmapisec.pl <-i bnx-file> <-a first-xmap-file> <-b 2nd-xmap-file>
# script version:1.0 (2016-06-21)
# Additional optional parameters are:
# -n <prefix for the output files> (default='isec_')>
# -c <minimal confidence score to be considered (default='undef')>
# -z zip the output to save space> (default OFF)>
# <-h to display this help>
```

### **mapisec.pl**
*[[General-tools](#general-tools)]*  

**REM: This script was first developped but should be preferred the upper one using xmap data as input because '.map' files are doomed to disappear in a future releases**

The perl script **[mapisec.pl](general-tools/mapisec.pl)** takes information from two MQR runs (from the 'MoleculeQualityReport.map' file) to split the corresponding BNX file into molecules that align to either, both, or none of the reference cmaps used for eather MQR. This allows creating subset of a BNX file that may be more specific for one or another reference assembly (which could for example represent parental genomes for a diploid) and make the resulting BNX data accessible for other applications like denovo assembly. An optional parameter allows filtering alignments by their 'Alignment Score' (warning, this score can take negative values!) to create more stringent datasets. The MQR runs may be performed using a lower than normal '-T' value in order to keep molecules that align with some degree of divergence (parental genomes are not necessariy identical to the haploid compound of a hybrid genome). *We would welcome your feedback after using this tool and reports of success would be a great reward for the work put into this script.*
```
Aim: Identify molecules specific to two ref-cmaps, ubiquitous, or not-aligning

# Usage: mapisec.pl <-i bnx-file> <-a first-map-file> <-b 2nd-map-file>
# script version:1.0.1 (2016-06-21)
# Additional optional parameters are:
# -n <prefix for the output files> (default='isec_')>
# -c <minimal confidence score to be considered (default='undef')>
# -z zip the output to save space> (default OFF)>
# <-h to display this help>
```

### **bnx0convert.pl**
*[[General-tools](#general-tools)]*  

The perl script **[bnx0convert.pl](general-tools/bnx0convert.pl)** converts old BNX version 0.1 data to teh curren tversion 1.2 format. It adds requierd fields with arbitrary values and brings your old data to compatibility with the current IrysView toolshed. Only one label is supported in this script as in all other scripts presented here. If you need to handle more than just one label, you will need to adapt the code by yourself.

```
Aim: Reformat old BNX 0.1 format to current version 1.2. Arbitrary values are used for backbone SNR and average intensity. You must provide a BNX file with -i
# script version:1.0 (2016-06-13)
## Usage: bnx0convert.pl <-i bnx-file>
# <-o run name (default 'UnknownRun'>
# <-h to display this help>
````

### **mqr2bnx.pl**
*[[General-tools](#general-tools)]*  

The perl script **[mqr2bnx.pl](general-tools/mqr2bnx.pl)** uses the xmap and input BNX from a MQR (quasi reference alignment obtained from BNG **RefAligner**) to identify BNX records that show homology to a reference (or genomic locus) and extract them to new BNX file. A minimal mapping confidence can be set to obtain BNX data of higher confidence. Finally, non-alignining molecules can also be saved to a second BNX file. The resulting BNX's can be denovo assembled or used as you wish them to be. Playing with MQR settings should allow producing datasets of variable specificity (to be tested :-)
```bash
# script version:1.01 (06-2016)
# Usage: mqr2bnx.pl <-b bnx-file> <-x xmap-file>
# Optional parameters:
# -c <minimal confidence score (default=0)>
# -n <save non-aligning BNX records to a second file (default OFF)>
# <-h to display this help>
```

### **bnxclnheader.pl**
*[[General-tools](#general-tools)]*  

The perl script **[bnxclnheader.pl](general-tools/bnxclnheader.pl)** will clean/shorten file path in the '# Run Data' lines of a BNX header.
```bash
## Usage: bnxclnheader.pl <-i bnx-file>
# Additional optional parameters are:
# <-z zip results (default OFF)>
# <-h to display this help>
```

### **bnxreheader.pl**
*[[General-tools](#general-tools)]*  

The perl script **[bnxreheader.pl](general-tools/bnxreheader.pl)** will swap the BNX header of a badly formatted BNX file with a correct header from a related BNX file. This script was made to correct multiple syntax issues in headers originated from IrysView 2.4 filtering or masking of repeats.
```bash
## Usage: bnxreheader.pl <-i bnx-file> <-t template-file>
# Additional optional parameters are:
# <-z zip results (default OFF)>
# <-h to display this help>
```

### **bedRename.pl**
*[[General-tools](#general-tools)]*  

The perl script **[bedRename.pl](general-tools/bedRename.pl)** will create a new BED file from a public file and replace the original chromosome names with the BNG translation provided with a key file (first column='official-name', second column='BNG-key' from Knicker). The resulting file can be viewed in **[IGV](https://www.broadinstitute.org/igv/)** together with BNG data.
```bash
## Usage: 
bedRename.pl <-i bed file (required)> <-k key file (required)>
# <-h to display this help>
```

### **cmap2bed.pl**
*[[General-tools](#general-tools)]*  

The perl script **[cmap2bed.pl](general-tools/cmap2bed.pl)** will create a BED file from a data.cmap file. The resulting file can be used with **[BEDTools](http://bedtools.readthedocs.org/en/latest/)** to go further.
```bash
Aim: Convert cmap data to BED5. You must provide a cmap file with -i
# Usage: cmap2bed.pl <-i cmap-file> 
# Optional parameters:
# -s <field number for BED-score (0-based; default to coverage=7)>
# (v1.2) 1:ContigLength 2:NumSites 3:SiteID 4:LabelChannel 
#        5:Position
#	 6:StdDev 7:Coverage 8:Occurrence 
#        9:GmeanSNR 10:lnSNRsd
# <-h to display this help>
```

### **xmap2bed.pl**
*[[General-tools](#general-tools)]*  

The perl script **[xmap2bed.pl](general-tools/xmap2bed.pl)** will create a BED5 file from a 'xmap' file. Only the part of the query that aligns to the reference is extracted and the result is expressed in REF(=anchor)-coordinates. Note that a size difference between query and reference matching regions will not be represented since reference coordinates cannot be modified. The resulting file can be viewed in **[IGV](https://www.broadinstitute.org/igv/)** or used with **[BEDTools](http://bedtools.readthedocs.org/en/latest/)** to go further. Users can filter and keep only alignments with a confidence greater than a given threshold (-x).
```bash
Aim: Convert xmap data to BED5. You must provide a xmap file with -i
# Usage: xmap2bed.pl <-i xmap-file>
# Optional parameters (v0.2) :
# -x <minimal value for score (default=0)>
# -c <coordinate system used <'q'=query/'r'=ref> (default='r')
# -n <field number for BED-name (1-based; default to SmapEntryID=1)>
#        1:XmapEntryID 2:QryContigID 3:RefcontigID1 4:QryStartPos 5:QryEndPos
#        6:RefStartPos 7:RefEndPos 8:Orientation 9:Confidence 10:HitEnumType 
#		11:Qrylen 12:Reflen 13:LabelChannel 14:Alignment
# -s <field number for BED-name (1-based; default to SmapEntryID=1)>
#        1:XmapEntryID 2:QryContigID 3:RefcontigID1 4:QryStartPos 5:QryEndPos
#        6:RefStartPos 7:RefEndPos 8:Orientation 9:Confidence 10:HitEnumType 
#		11:Qrylen 12:Reflen 13:LabelChannel 14:Alignment
# -k <key file (when provided, will rename the sequences to their original naming (default absent)>
# <-h to display this help>
```

### **xmap2bed12.pl**
*[[General-tools](#general-tools)]*  

The perl script **[xmap2bed12.pl](general-tools/xmap2bed12.pl)** will create a BED12 file from a 'xmap' file. The part of the query that aligns to the reference is represented in thick block while additional query ends not matching the reference are produced as thin blocks left and right from the match. Note that a size difference between query and reference matching regions will not be represented since reference coordinates cannot be modified. The resulting file can be viewed in **[IGV](https://www.broadinstitute.org/igv/)** or used with **[BEDTools](http://bedtools.readthedocs.org/en/latest/)** to go further. Users can filter and keep only alignments with a confidence greater than a given threshold (-x). See some example data **[here](xmap2bed12.pl_example.rmd)**.
```bash
Aim: Convert xmap data to BED12. You must provide a xmap file with -i
# script version:1.1 (05-2016)
# Usage: xmap2bed12.pl <-i xmap-file>
# Optional parameters (xmap v0.2) :
# -x <minimal value for score (default=0)>
# -r <RGB feature color 255,0,0=red (default=0 | black)>
# -k <key file (when provided, will rename the sequences to their original naming (default absent)>
# -v <report verbose summary>
# <-h to display this help>
```

### **smap2bed.pl**
*[[General-tools](#general-tools)]*  

The perl script **[smap2bed.pl](general-tools/smap2bed.pl)** will create a BED file from a 'smap' file. The resulting file can be used with **[IGV](https://www.broadinstitute.org/igv/)** or **[BEDTools](http://bedtools.readthedocs.org/en/latest/)** to go further. Note that we need to work more on that piece of code as the SV data is not straightforward. Any suggestions would be welcome to make this one more useful. The script also reports the distribution of Confidence scores and allows identifying relevant cutoff values for a second run (it is possible to get the cutoff value for any given percentile limit by modifying the '-p' argument).  See some example data **[here](smap2bed.pl_example.rmd)**.
```bash
# Usage: smap2bed.pl <-i smap-file>
# Optional parameters (smap v0.4) :
# -x <minimal value for Confidence score (default=0, '-1' is used for complex calls)>
# -c <coordinate system used <'q'=query/'r'=ref> (default='r')
# -n <field number for BED-name (1-based; default to SmapEntryID=1)>
#        1:SmapEntryID 2:QryContigID 3:RefcontigID1 4:RefcontigID2 5:QryStartPos 6:QryEndPos
#        7:RefStartPos 8:RefEndPos 9:Confidence 10:Type 11:XmapID1 12:XmapID2 13:LinkID
#       14:QryStartIdx 15:QryEndIdx 16:RefStartIdx 17:RefEndIdx
# -s <field number for BED-score (1-based; default to Confidence=9)>
#        1:SmapEntryID 2:QryContigID 3:RefcontigID1 4:RefcontigID2 5:QryStartPos 6:QryEndPos
#        7:RefStartPos 8:RefEndPos 9:Confidence 10:Type 11:XmapID1 12:XmapID2 13:LinkID
#       14:QryStartIdx 15:QryEndIdx 16:RefStartIdx 17:RefEndIdx
# -p <percentile for Confidence distribution (default=95>)
# -k <key file (when provided, will rename the sequences to their original naming (default absent)>
# <-h to display this help>
```

### **bnxsplitter.pl**
*[[General-tools](#general-tools)]*  

The perl script **[bnxsplitter.pl](general-tools/bnxsplitter.pl)** will split data from a BNX file (or archive thereof) into five separate 'invisible' TSV files available for down-processing using **R** (or your favorite script).

* .header.tsv
* .zero.tsv
* .one.tsv
* .qx11.tsv
* .qx12.tsv

The files are created next to the input file and are made invisible with a starting -dot- This can easily be changed in the code.
```bash
Aim: Split a BNX file into its components. You must provide a BNX file with -i
## Usage: bnxsplitter.pl <-i bnx-file>
# <-h to display this help>
```

### **labeldensity.pl**
*[[General-tools](#general-tools)]*  

The perl tool **[labeldensity.pl](general-tools/labeldensity.pl)** will find regions of the genome that show abnormal label densities (none or high). 

To achieve this, it proceeds as follows:

* Search for selected nicking enzyme sites in a multifasta and save results in BED format (depends on: **[restrict2bed.pl](https://github.com/BITS-VIB/ngs-tools/blob/master/fasta-tools/restrict2bed.pl)**)
* Create the chromosome / contig-size list in BED format a from multifasta (depends on: **[fasta2chromsizes.pl](https://github.com/BITS-VIB/ngs-tools/blob/master/fasta-tools/fasta2chromsizes.pl)**)
* Create windows for each chromosome / contig in teh former file and save in BED format (required: **bedtools makewindows**)
* Compare the first and last BED files and record the count of nicking sites in each bin (required: **bedtools map**)
* Sort the obtained BED file naturally (requires a recent version of **GNU sort**) for visualisation (eg. **[IGV](https://www.broadinstitute.org/igv/)**).

```bash
## Usage: labeldensity.pl <-i fasta-file> <-n 'nicker(s)'>
# multiple allowed separated by ',')>
#  'Nt-BspQI' => 'GCTCTTC',
#  'Nt-BbvCI' => 'CCTCAGC',
#  'Nb-BsMI'  => 'GAATGC',
#  'Nb-BsrDI' => 'GCAATG'
# Additional optional parameters are:
# <-t title ('label-density')>
# <-l minimal length for dna sequence (20000)>
# <-b bin width for computing label density (100000)>
# <-h to display this help>
```

### **cmap2renum.pl**
*[[General-tools](#general-tools)]*  

The perl tool **[cmap2renum.pl](general-tools/cmap2renum.pl)** takes one reference cmap and its matching key-file generated by 'Knicker' and renumbers all cmaps starting from 1 in both files. A new pair of files is saved to disk with added prefix. Such operation is required when the original cmap contained high values for the cmap IDs (over 100,000) which is not supported by downstream steps like hybrid scaffolding). large ID numbers may come from very large contig lists where a number of sequences have been filtered out due to Knicker cutoffs, leaving holes in the ID range and breaching the limit of 100,000.

```bash
# Usage: cmap2renum.pl <-c cmap-file> <-k key-file>
# Optional parameters:
# -p <prefix ('new_')>
# -h <this help message>
```

In order to clean your assembly file, you may consider applying the next two perl scripts before using 'Knicker'.

### **findNregions.pl**
*[[General-tools](#general-tools)]*  

The perl tool **[findNregions.pl](general-tools/findNregions.pl)** find regions of N's from a reference multi-fasta file and the corresponding knicker key table. It stores the coordinate of all hits to BED for loading in IrysView as track with sequence titles renamed using the key file. Such track may prove useful to identify issues associated with sequence gaps of incorrect size introduced in assemblies.

```bash
## Usage: findNregions.pl <-i fasta-file> <-k key-file to rename contigs>
# Additional optional parameters are:
# <-l minsize in bps (100)>
# <-h to display this help>
```

### **fastaFiltLength.pl**
*[[General-tools](#general-tools)]*  

The BIO-perl script **[fastaFiltLength.pl](fasta-tools/fastaFiltLength.pl)** will filter a multifasta file and keep only sequence with length &gt; min and $lt; max values. Was created to filter genome assemblies containing multiple small files.

```bash
## Usage: fastaFiltLength.pl <-i fasta_file (required)>
# Additional optional parameters are:
# <-o outfile_name (filtered_)>
# <-m minsize (undef)>
# <-x maxsize (undef)>
# <-h to display this help>
```

### **fastaSortLength.pl**
*[[General-tools](#general-tools)]*  

The BIO-perl script **[fastaSortLength.pl](fasta-tools/fastaSortLength.pl)** will sorts a multifasta file by decreasing or increasing order. Was created to clean input fasta files before applying Knicker (BionanoGenomics).
```bash
## Usage: fastaSortlength.pl <-i fasta-file> <-o size-order ('i'=increasing | 'd'=decreasing)>
# <-h to display this help>
```

### **fastaRename.pl**
*[[General-tools](#general-tools)]*  

The BIO-perl script **[fastaRename.pl](fasta-tools/fastaRename.pl)** will rename headers of a multifasta file using the index file generated by Knicker (BionanoGenomics). The resulting file can be used with data exported from BNG in **[IGV](https://www.broadinstitute.org/igv/)**.
```bash
## Usage: 
fastaRename.pl <-i fasta_file (required)> <-k key file (required)>
# <-h to display this help>
```

### **fastaPad20k.pl**
*[[General-tools](#general-tools)]*  

The perl script **[fastaPad20k.pl](general-tools/fastaPad20k.pl)** adds N's and a nicking motif to sequences shorter than 20k bases. Not sure it will rescue any but this was the idea.

```bash
## Usage: fastaPad20k.pl
## script version:1.0
# <-i fasta_file (required)>
# <-n nicking motif (required)>
#   Nt.BspQI: GCTCTTCN
#   Nt.BbvCI: CCTCAGC
#   Nb.Bsml: GAATGCN
#   Nb.BbvCI: CCTCAGC
#   Nb.BsrD1: GCAATGNN
#   Nb.BssSI: CACGAG
# <-o outfile | default from infile name>
# <-z zip results (default OFF)>
# <-h to display this help>
```

### **FastaExtractList.pl**
*[[General-tools](#general-tools)]*  

The perl script **[FastaExtractList.pl](general-tools/FastaExtractList.pl)** extracts a list of Fasta records from a multifasta file into a new file.

```bash
## Usage: FastaExtractList.pl 
## script version:1.0
# <-i fasta_file (required)> 
# <-r record_list file (required)>
# <-o outfile | (optional) default from infile name>
# <-z zip results (default OFF)>
# <-h to display this help>
```

### **CmapExtractList.pl**
*[[General-tools](#general-tools)]*  

The perl script **[CmapExtractList.pl](general-tools/CmapExtractList.pl)** extracts a list of cmap records from a cmap file into a new file. Can be used to create subset of a optical maps from sequence names returned in a Hybrid-scaffold-not-scaffolded.fasta file.

```bash
## Usage: CmapExtractList.pl 
## script version:1.0
# <-i cmap_file (required)> 
# <-r record_list file (required)>
# <-k key file (required)>
# <-o outfile | (optional) default from infile name>
# <-h to display this help>
```

### **run_SV.sh**
*[[General-tools](#general-tools)]*  

The bash script **[run_SV.sh](general-tools/run_SV.sh)** runs SV analysis from files obtained in denovo assembly.

```bash
# Usage: run_SV.sh -r <path to the reference.cmap> -i <assembly-folder>
# script version 1.0, 2016_08_16
# [optional: -t <path to Refaligner (default to "/home/bionano/tools")>]
# [optional: -q <query_path containing assembly cmaps (default to <assembly>/output/contigs/exp_refineFinal1)>]
# [optional: -o <output folder (defaults in input map folder with same name + suffix "_sv")>]
# [optional: -p <path to the python script "runSV.py" (default to "/home/bionano/scripts")>]
# [optional: -a <optArgument.xml (default to "<assembly-folder>/optArguments_XXX.xml")>]
# [optional: -T <number of threads to use (default=4)>]
# [optional: -j <threads per job (default=4)>]
# [optional: -b <BED file with GAPs (unset by default)>]
# [optional: -e <.err file (defaults to <query_path>/alignref_final/EXP_REFINEFINAL1.err)>]
# [optional: -E <.errbin file (defaults to <query_path>/alignref_final/EXP_REFINEFINAL1.errbin)>]
# [optional: -C <CXML file for running on cluster (unset by default)>]
# [optional: -s <SV job configuration 0="single job" (default), 1="single job per contig" (not recommended), 2="grouped">
# [-h for this help]
```

### **run_HS.sh**
*[[General-tools](#general-tools)]*  

The bash script **[run_HS.sh](general-tools/run_HS.sh)** runs a hybridScaffold analysis from files obtained in denovo assembly. The manual hybridscaffold using edited conflicts (-M) is not yet supported by this code which does only the regular HS. Please check that the default used xml files match your genome properties or link to the correct files.

```bash
# Usage: run_HS.sh
# script version 1.2, 2017_01_15
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
# [-h for this help]
```

### **run_DNloc.sh**
*[[General-tools](#general-tools)]*  

The bash script **[run_DNloc.sh](general-tools/run_DNloc.sh)** runs a denovo assembly using only local cpu resources (for servers without Xeon-Phi cards). Please read the output of **python ${SCRIPTS}/pipelineCL.py -help**  for more arguments. 

REM: The denovo assembly is here done with a reference cmap and computing autonoise against it. When your reference is not that good, please run a manual command instead of this script and omit the reference (-r) and the autonoise (-y) arguments. You can also limit the refinements to -i 3 to obtain a rough unsupervised assembly that can be used as reference for autonoise in a subsequent denovo run.

```bash
# Usage: run_DNloc.sh
# script version 1.0, 2016_11_06
## arguments
# [required: -b <molecule BNX file to assemble>]
# [required: -r <BioNano ref CMAP file for noise computation and stats>]
# [required: -x <optArgument.xml>]
# [optional: -o <assembly-base-folder (default current folder)>]
# [optional: -s <pipelineCL.py path (required if not in the default location)]
# [optional: -t <max-threads | 8 >]
# [optional: -j <max-jobs (max-thread/2) | 4 >]
# [-h for this help]
```

<details>
<summary>**BioNano Genomics ManPage for pipelineCL.py** (*click to expand*)</summary>
```bash
usage: pipelineCL.py [-h] [-T T] [-j MAXTHREADS] [-N N] [-G BED] [-i ITER]
                     [-I IMG] [-b BNX] [-l LOCAL] [-t TOOLS] [-B BYPASS]
                     [-e EXP] [-r REF] [-s GROUPSV] [-n] [-x] [-c CLEANUP]
                     [-g SIZE] [-C CXML] [-w] [-a XML] [-L LAMBDAREF]
                     [-p PERF] [-d] [-f NGSARG] [-F] [-u] [-U [GROUPCONTIGS]]
                     [-v [VERSION]] [-S] [-V RUNSV] [-A] [-y] [-Y] [-m] [-H]

Pipeline for de novo assembly - BioNano Genomics

optional arguments:
  -h, --help         show this help message and exit
  -T T               Available threads per Node [default 1]
  -j MAXTHREADS      Threads per job [default 1]
  -N N               Number of split bnx files; number of pairwise jobs is
                     N*(N-1)/2 (optional, default: -T)
  -G BED             Bed file for gaps, used in structural variation (SV)
                     detection to check for SV overlap with reference gaps
  -i ITER            Number of extension and merge iterations (default=1, must
                     be in range [0,10], use 0 to skip)
  -I IMG             File with listed paths for image processing; no longer
                     supported--do not use (use without .bnx)
  -b BNX             Input molecule (.bnx) file, required
  -l LOCAL           Location of output files root directory, required, will
                     be created if does not exist; if does exist, will
                     overwrite contents (may be error-prone)
  -t TOOLS           Location of executable files (RefAligner and Assembler,
                     required)
  -B BYPASS          Skip steps, using previous result. <= 0:None,
                     1:ImgDetect, 2:NoiseChar/Subsample, 3:Pairwise,
                     4:Assembly, 5:RefineA, 6:RefineB, (7:RefineNGS -f only),
                     7:merge0, 8+(i-1)*2:Ext(i), 9+(i-1)*2:Mrg(i),
                     N+1:alignmol
  -e EXP             Output file prefix (optional, default = exp)
  -r REF             Reference file (must be .cmap), to compare resulting
                     contigs (optional)
  -s GROUPSV         SV jobs configuration: 0 = single job (required for
                     correct haplotype calls), 1 = single job per contig (not
                     recommended), 2 = grouped (default 0; optional)
  -n                 Evaluate single molecule noise characterization
  -x                 Exit after auto noise (noise characterization), do not
                     preform de novo assembly
  -c CLEANUP         Remove contig results (0 - keep all (default), 1 - remove
                     intermediate files, 2 - store in sqlite, 3 - store in
                     sqlite and remove)
  -g SIZE            Organism genome size estimate in megabases, used for
                     tuning assembly parameters [optional, if > 0, will modify
                     parameters, if == 0, ignored, must be float]
  -C CXML            Run on cluster, read XML file for submission arguments
                     (optional--will not use cluster submission if absent)
  -w                 Wipe clean previous contig results
  -a XML             Read XML file for parameters (required)
  -L LAMBDAREF       Lambda phage is spiked in, used for molecule scaling
                     (only used with -I input)
  -p PERF            Log performance in pipelineReport 0=None, 1=time, 2=perf,
                     3=time&perf (default=1)
  -d                 Retired option (contig subdirectories), always enabled.
  -f NGSARG          Directory which contains NGS contigs as cmaps for use in
                     refineNGS (no longer supported).
  -F                 When using -f, disable min contigs check for all stages
                     prior to refineNGS (no longer supported).
  -u                 Do not perform final refinement (not recommended).
  -U [GROUPCONTIGS]  Group contigs in refinement and extension stages [default
                     ON, use 0 to disable]
  -v [VERSION]       Print version; exit if argument > 1 supplied.
  -S                 Log stdout/stderr of all RefAligner calls in files with
                     same names as output files and suffix .stdout (default on
                     --use this flag to turn off).
  -V RUNSV           Detect structural variations. Default: only after final
                     stage (normally refineFinal); if argument 2, also after
                     refineB; if argument 0, disable.
  -A                 Align molcules to final contigs (ON by default, use this
                     to turn off).
  -y                 Automatically determine noise parameters (requires
                     reference; optional, default off)
  -Y                 Disable scan scaling in auto noise (default on with
                     reference and -y)
  -m                 Disable molecule vs reference alignments (default on with
                     reference)
  -H                 Use HG19 (human genome) as reference, loaded from
                     Analysis/SV/CopyNumberProfiles/DATA. Overrides -r
                     argument.
```
</details>

## SysAdmin-tools
*[[back-to-top](#top)]*  

Those additional tools that we had to develop to troublechoot problems.

### **logphicards.sh**
*[[SysAdmin-tools](#sysadmin-tools)]*  

The perl script **[logphicards.sh](sysadmin/logphicards.sh)** logs several metrics for the 6 Xeon-phi cards present in our server and stores the results in a text file. The log file is then used to plot the different parameters as shown in a demo report attached **[here](sysadmin/thinkmate_logging.pdf)**.
```bash
# Usage: logphicards.sh
#    -t <log-frequency in sec (default 60sec)>
```

## Shiny-apps
*[[back-to-top](#top)]*  

Those additional tools belong on a Shiny server and will execute R code in a interactive manner (please refder to https://shiny.rstudio.com/ for info about Shiny).

### **BNX_viewer.shinyapp** 
*[[Shiny-apps](#shiny-apps)]*

The **[BNX_viewer.shinyapp](shiny-apps/BNX_viewer.shinyapp)** app loads BNX data and creates scattreplots with filtering (no export at this point, one can easily use built-in functions for that purpose). A live version was posted to https://nucleomics-core.shinyapps.io/bnx_filter/ (this version is size limited, please use only small BNX files to keep this server alive). A sample BNX file is present in the 'Data' subfolder for your convenience [(sample.bnx)](https://github.com/Nucleomics-VIB/bionano-tools/raw/master/shiny-apps/BNX_viewer.shinyapp/Data/sample.bnx).


*[[back-to-top](#top)]*  

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).
