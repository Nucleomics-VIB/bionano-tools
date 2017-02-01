#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Dumper;

# parse assembly report 'exp_informaticsReport.txt'
# reformat and output nicer 1-page
#
# Stephane Plaisance (VIB-NC+BITS) 2016/12/08; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

getopts('i:h');
our($opt_i, $opt_h);

my $usage="## Usage: assemblyReport.pl <-i exp_informaticsReport.txt>
# <-h to display this help>";

####################
# declare variables
####################

my $report = $opt_i || die $usage."\n";
defined($opt_h) && die $usage."\n";

# slurp data in variable
open (DATA, $report);
my $input;
my @content;
while (<DATA>){
       $input .= $_;
       }
# convert to array of text blocks
@content = split(/\n{2,}/, $input); # note the new pattern

# parse blocks by type
my $command = shift @content;
## split command arguments for better readability
my @cmdsplit = split(" -", $command);
my $cmdout = join("\n\t-", @cmdsplit);
my $title = shift @content;

# start output
print STDOUT $title."\n-----------\n#Command :".$cmdout."\n";

# handle other text blocks
my $statreg = 'Reading molecule stats from';
my $noisereg = 'noise parameters';
my $molecreg = 'Molecules Aligned';
my $stagereg = 'Stage Summary:';
my $svreg = 'SV detect:';

# counter and array for stages
my $id=0;
my @stagedata = ();

foreach (@content) {
	if ($_ =~ m/\Q$statreg/)
		{
		my @rows = split("\n", $_);
		my $infile = shift @rows;
		$infile =~ s/\Q$statreg//;
		$infile =~ s/\:$//;
		my $title = shift @rows;
		print STDOUT "--------\n".$title.$infile."\n";
		foreach (@rows) {
			my @fields = split(": ", $_);
			map {s/\ +/\ /g; } @fields;
			print STDOUT join("\t", @fields)."\n";
			}
		}

	elsif ($_ =~ m/\Q$noisereg/)
		{
		my @rows = split("\n", $_);
		my $title = shift @rows;
		print STDOUT "--------\n".$title."\n";
		foreach (@rows) {
			my @fields = split(":", $_);
			map {s/\ +/\ /g; } @fields;
			print STDOUT join("\t", @fields)."\n";
			}
		}

	elsif ($_ =~ m/\Q$molecreg/)
		{
		my @rows = split("\n", $_);
		my $title = shift @rows;
		print STDOUT "--------\n".$title."\n";
		foreach (@rows) {
			my @fields = split(": ", $_);
			map {s/\ +/\ /g; } @fields;
			print STDOUT join("\t", @fields)."\n";
			}
		}

	elsif ($_ =~ m/\Q$stagereg/)
		{
		# new stage block
		$id++;

		my @rows = split("\n", $_);
		my $title = shift @rows;
		my $label = $title =~ s/\QStage Summary: //r;

		# name the row in stagedata
		$stagedata[$id]{'stage'} = $label;

		print STDOUT "--------\n".$title."\n";

		foreach (@rows) {
			my @fields = split(": ", $_);
			map {s/\ +/\ /g; } @fields;
			print STDOUT join("\t", @fields)."\n";
			$stagedata[$id]{$fields[0]} = $fields[1];
			}
		}

	elsif ($_ =~ m/\Q$svreg/)
		{
		my @rows = split("\n", $_);
		my $infile = shift @rows;
		my $trash = shift @rows;
		print STDOUT "--------\n".$infile."\n";
		my @table = grep { $_ !~ /^N / } @rows;
		my @tot = grep { /^N / } @rows;
		# print table
		foreach (@table) {
			my @fields = split(": ", $_);
			map {s/\ +/\ /g; } @fields;
			map {s/^\ +//g; } @fields;
			print STDOUT join("\t", @fields)."\n";
			}
		# print totals
		print STDOUT "# total counts:\n";
		foreach (@tot) {
			my @fields = split(": ", $_);
			print STDOUT join("\t", @fields)."\n";
			}
		}

	else
		{
		print "## unexpected block\n";
		}
	}

exit 0;

#########################
#### example content ####

# /home/bionano/scripts//pipelineCL.py -U -d -T 228 -j 228 -N 4 -i 5 -a /home/bionano/data/test_new/2258lizard-FM100k/optArguments_haplotype.xml -w -t /home/bionano/tools/ -l /home/bionano/data/test_new/2258lizard-FM100k/output -b /home/bionano/data/test_new/2258lizard-FM100k/Molecules.bnx -r /home/bionano/data/test_new/2258lizard-FM100k/HLtupMer3_BspQI.cmap -y -C /home/bionano/scripts//clusterArguments.xml
#
# Informatics Report Version: 2.3.0
#
# Reading molecule stats from /home/bionano/data/test_new/2258lizard-FM100k/output/all.bnx:
# Molecule Stats:
# N mols: 1539583
# Total len (Mb): 264672.976
# Avg len (kb)  :    171.912
# Mol N50 (kb)  :    174.448
# Lab (/100kb)  :     10.813
#
# Automatically determined noise parameters:
# FP:0.968163
# FN:0.152197
# sf:0.212651
# sd:-0.0985
# sr:0.04102
# bpp:511.02
# readparameters:/home/bionano/data/test_new/2258lizard-FM100k/output/contigs/auto_noise/autoNoise1.errbin
#
# Stage Summary: CharacterizeDefault Assembly
# N Genome Maps: 4253
# Total Genome Map Len  (Mb): 2556.890
# Avg. Genome Map Len   (Mb):    0.601
# Median Genome Map Len (Mb):    0.415
# Genome Map n50        (Mb):    0.813
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.284
# N Genome Maps total align      : 3366 (0.79)
# Total Aligned Len (Mb)            : 2067.952
# Total Aligned Len / Ref Len       :    1.038
# Total Unique Aligned Len (Mb)     : 1793.699
# Total Unique Aligned Len / Ref Len:    0.900
#
# Reading molecule stats from /home/bionano/data/test_new/2258lizard-FM100k/output/contigs/auto_noise/autoNoise1_rescaled.bnx:
# Molecule Stats:
# N mols: 775015
# Total len (Mb): 169363.861
# Avg len (kb)  :    218.530
# Mol N50 (kb)  :    213.551
# Lab (/100kb)  :     10.516
# Ref    Cov (x):     85.026
#
# Molecules Aligned to Reference:
# N mol align       :    449108
# Mol fraction align:         0.579
# Tot align len (Mb):     85845.0
# Effective Cov (x) :        43.527
# Avg align len (kb):       191.1
# Fraction align len:         0.507
# Tot confidence    :   7439689.5
# Avg confidence    :        16.6
# Avg FP(/100kb)    :         0.98
# Avg FP ratio      :         0.070
# Avg FN ratio      :         0.166
# Avg bpp           :       511.1
# Avg sf            :         0.221
# Avg sd            :        -0.103
# Avg sr            :         0.042
#
# Stage Summary: CharacterizeDefault refineB1
# N Genome Maps: 4422
# Total Genome Map Len  (Mb): 2326.832
# Avg. Genome Map Len   (Mb):    0.526
# Median Genome Map Len (Mb):    0.344
# Genome Map n50        (Mb):    0.772
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.168
# N Genome Maps total align      : 3716 (0.84)
# Total Aligned Len (Mb)            : 2063.146
# Total Aligned Len / Ref Len       :    1.036
# Total Unique Aligned Len (Mb)     : 1831.913
# Total Unique Aligned Len / Ref Len:    0.920
#
# Stage Summary: CharacterizeDefault Merge0
# N Genome Maps: 3290
# Total Genome Map Len  (Mb): 2065.055
# Avg. Genome Map Len   (Mb):    0.628
# Median Genome Map Len (Mb):    0.463
# Genome Map n50        (Mb):    0.856
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.037
# N Genome Maps total align      : 2813 (0.86)
# Total Aligned Len (Mb)            : 1858.813
# Total Aligned Len / Ref Len       :    0.933
# Total Unique Aligned Len (Mb)     : 1830.592
# Total Unique Aligned Len / Ref Len:    0.919
#
# Stage Summary: CharacterizeDefault extension1_1
# N Genome Maps: 3119
# Total Genome Map Len  (Mb): 2034.846
# Avg. Genome Map Len   (Mb):    0.652
# Median Genome Map Len (Mb):    0.486
# Genome Map n50        (Mb):    0.868
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.022
# N Genome Maps total align      : 2779 (0.89)
# Total Aligned Len (Mb)            : 1851.068
# Total Aligned Len / Ref Len       :    0.929
# Total Unique Aligned Len (Mb)     : 1828.112
# Total Unique Aligned Len / Ref Len:    0.918
#
# Stage Summary: CharacterizeDefault Merge1
# N Genome Maps: 3061
# Total Genome Map Len  (Mb): 2022.861
# Avg. Genome Map Len   (Mb):    0.661
# Median Genome Map Len (Mb):    0.496
# Genome Map n50        (Mb):    0.873
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.016
# N Genome Maps total align      : 2744 (0.90)
# Total Aligned Len (Mb)            : 1844.214
# Total Aligned Len / Ref Len       :    0.926
# Total Unique Aligned Len (Mb)     : 1828.036
# Total Unique Aligned Len / Ref Len:    0.918
#
# Stage Summary: CharacterizeDefault extension1_2
# N Genome Maps: 3019
# Total Genome Map Len  (Mb): 2014.032
# Avg. Genome Map Len   (Mb):    0.667
# Median Genome Map Len (Mb):    0.499
# Genome Map n50        (Mb):    0.874
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.011
# N Genome Maps total align      : 2734 (0.91)
# Total Aligned Len (Mb)            : 1842.248
# Total Aligned Len / Ref Len       :    0.925
# Total Unique Aligned Len (Mb)     : 1827.525
# Total Unique Aligned Len / Ref Len:    0.917
#
# Stage Summary: CharacterizeDefault Merge2
# N Genome Maps: 3009
# Total Genome Map Len  (Mb): 2012.136
# Avg. Genome Map Len   (Mb):    0.669
# Median Genome Map Len (Mb):    0.500
# Genome Map n50        (Mb):    0.875
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.010
# N Genome Maps total align      : 2731 (0.91)
# Total Aligned Len (Mb)            : 1841.925
# Total Aligned Len / Ref Len       :    0.925
# Total Unique Aligned Len (Mb)     : 1827.555
# Total Unique Aligned Len / Ref Len:    0.917
#
# Stage Summary: CharacterizeDefault extension1_3
# N Genome Maps: 3001
# Total Genome Map Len  (Mb): 2009.582
# Avg. Genome Map Len   (Mb):    0.670
# Median Genome Map Len (Mb):    0.501
# Genome Map n50        (Mb):    0.875
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.009
# N Genome Maps total align      : 2728 (0.91)
# Total Aligned Len (Mb)            : 1841.102
# Total Aligned Len / Ref Len       :    0.924
# Total Unique Aligned Len (Mb)     : 1827.039
# Total Unique Aligned Len / Ref Len:    0.917
#
# Stage Summary: CharacterizeDefault Merge3
# N Genome Maps: 2998
# Total Genome Map Len  (Mb): 2008.934
# Avg. Genome Map Len   (Mb):    0.670
# Median Genome Map Len (Mb):    0.501
# Genome Map n50        (Mb):    0.876
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.009
# N Genome Maps total align      : 2727 (0.91)
# Total Aligned Len (Mb)            : 1840.925
# Total Aligned Len / Ref Len       :    0.924
# Total Unique Aligned Len (Mb)     : 1827.039
# Total Unique Aligned Len / Ref Len:    0.917
#
# Stage Summary: CharacterizeDefault extension1_4
# N Genome Maps: 2997
# Total Genome Map Len  (Mb): 2007.955
# Avg. Genome Map Len   (Mb):    0.670
# Median Genome Map Len (Mb):    0.501
# Genome Map n50        (Mb):    0.876
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.008
# N Genome Maps total align      : 2727 (0.91)
# Total Aligned Len (Mb)            : 1840.551
# Total Aligned Len / Ref Len       :    0.924
# Total Unique Aligned Len (Mb)     : 1826.693
# Total Unique Aligned Len / Ref Len:    0.917
#
# Stage Summary: CharacterizeDefault Merge4
# N Genome Maps: 2997
# Total Genome Map Len  (Mb): 2007.955
# Avg. Genome Map Len   (Mb):    0.670
# Median Genome Map Len (Mb):    0.501
# Genome Map n50        (Mb):    0.876
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.008
# N Genome Maps total align      : 2727 (0.91)
# Total Aligned Len (Mb)            : 1840.551
# Total Aligned Len / Ref Len       :    0.924
# Total Unique Aligned Len (Mb)     : 1826.693
# Total Unique Aligned Len / Ref Len:    0.917
#
# Stage Summary: CharacterizeDefault extension1_5
# N Genome Maps: 2996
# Total Genome Map Len  (Mb): 2007.157
# Avg. Genome Map Len   (Mb):    0.670
# Median Genome Map Len (Mb):    0.502
# Genome Map n50        (Mb):    0.875
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.008
# N Genome Maps total align      : 2727 (0.91)
# Total Aligned Len (Mb)            : 1840.388
# Total Aligned Len / Ref Len       :    0.924
# Total Unique Aligned Len (Mb)     : 1826.647
# Total Unique Aligned Len / Ref Len:    0.917
#
# Stage Summary: CharacterizeDefault Merge5
# N Genome Maps: 2996
# Total Genome Map Len  (Mb): 2007.157
# Avg. Genome Map Len   (Mb):    0.670
# Median Genome Map Len (Mb):    0.502
# Genome Map n50        (Mb):    0.875
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.008
# N Genome Maps total align      : 2727 (0.91)
# Total Aligned Len (Mb)            : 1840.388
# Total Aligned Len / Ref Len       :    0.924
# Total Unique Aligned Len (Mb)     : 1826.647
# Total Unique Aligned Len / Ref Len:    0.917
#
# Stage Summary: CharacterizeFinal refineFinal1
# Diploid N Genome Maps: 5552
# Diploid Genome Map Len        (Mb): 3809.186
# Diploid Avg. Genome Map Len   (Mb):    0.686
# Diploid Median Genome Map Len (Mb):    0.520
# Diploid Genome Map n50        (Mb):    0.896
# Haploid N Genome Maps: 2989
# Haploid Genome Map Len        (Mb): 1964.958
# Haploid Avg. Genome Map Len   (Mb):    0.657
# Haploid Median Genome Map Len (Mb):    0.491
# Haploid Genome Map n50        (Mb):    0.874
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.912
# N Genome Maps total align      : 5130 (0.92)
# Total Aligned Len (Mb)            : 3468.285
# Total Aligned Len / Ref Len       :    1.741
# Total Unique Aligned Len (Mb)     : 1796.706
# Total Unique Aligned Len / Ref Len:    0.902
#
# Stage Summary: CharacterizeDefault refineFinal1
# Diploid N Genome Maps: 5552
# Diploid Genome Map Len        (Mb): 3809.186
# Diploid Avg. Genome Map Len   (Mb):    0.686
# Diploid Median Genome Map Len (Mb):    0.520
# Diploid Genome Map n50        (Mb):    0.896
# Haploid N Genome Maps: 2989
# Haploid Genome Map Len        (Mb): 1964.958
# Haploid Avg. Genome Map Len   (Mb):    0.657
# Haploid Median Genome Map Len (Mb):    0.491
# Haploid Genome Map n50        (Mb):    0.874
# Total Ref Len   (Mb): 1991.901
# Total Genome Map Len / Ref Len :    1.912
# N Genome Maps total align      : 5131 (0.92)
# Total Aligned Len (Mb)            : 3519.240
# Total Aligned Len / Ref Len       :    1.767
# Total Unique Aligned Len (Mb)     : 1821.791
# Total Unique Aligned Len / Ref Len:    0.915
#
# SV detect: svdetect_exp_refineFinal1_sv
# map filename : len (Mb) : align len (Mb) (ratio) : N SV
#                  SV type :   N : N/tot
#                  complex :   36 : 0.0029
#                 deletion : 5135 : 0.4132
#            deletion_tiny :    6 : 0.0005
#                      end :  409 : 0.0329
#                insertion : 6677 : 0.5373
#           insertion_tiny :    2 : 0.0002
#                inversion :    7 : 0.0006
#        inversion_partial :    7 : 0.0006
#   translocation_interchr :  149 : 0.0120
# N contigs  :      1
# N sv total :  12428  (12428.0000 per contig)
# N sv filter:  11975  (11975.0000 per contig)
#
# Reading molecule stats from /home/bionano/data/test_new/2258lizard-FM100k/output/contigs/auto_noise/autoNoise1_rescaled.bnx:
# Molecule Stats:
# N mols: 775015
# Total len (Mb): 169363.861
# Avg len (kb)  :    218.530
# Mol N50 (kb)  :    213.551
# Lab (/100kb)  :     10.516
# Contig Cov (x):     44.462
#
# Molecules Aligned to Assembly:
# N mol align       :    523489
# Mol fraction align:         0.675
# Tot align len (Mb):    103776.5
# Effective Cov (x) :        27.708
# Avg align len (kb):       198.2
# Fraction align len:         0.613
# Tot confidence    :  11348413.5
# Avg confidence    :        21.7
# Avg FP(/100kb)    :         0.55
# Avg FP ratio      :         0.045
# Avg FN ratio      :         0.121
# Avg bpp           :       510.9
# Avg sf            :         0.174
# Avg sd            :        -0.061
# Avg sr            :         0.028

