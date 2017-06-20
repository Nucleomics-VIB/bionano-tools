#!/usr/bin/perl

# Print pretty TEXT report from BioNano 'exp_informaticsReport.txt'
# Takes haploid reports as well as diploid ones
# Stephane Plaisance (VIB-NC+BITS) 2017/06/20; v1.0
# visit our Git: https://github.com/BITS-VIB

use warnings;
use strict;
use Getopt::Std;
use File::Basename;

# specials
#use Data::Dumper;
use Array::Transpose;
# @array=transpose(\@array);

my $usage="## Usage: 
prettyReport.pl <-i exp_informaticsReport.txt> <-o base-name for output (optional)>
# <-h to display this help>";

####################
# declare variables
####################
getopts('i:o:h');
our ($opt_i, $opt_o, $opt_h);

my $infile = $opt_i || die $usage."\n";

my $outbase = $opt_o || "prettyReport";
defined($opt_h) && die $usage."\n";

open FILE, $infile or die $!;
my $outpath = dirname($infile);
my $outfile = $outpath."/".$outbase.".txt";

# declare variables
my $command;
my $pipelineversion;
my @stats;
my @statst;
my @stages;
my @stagest;
my @final;
my @finalt;
my @aligns;
my @alignst;
my @svres;

while (my $line = <FILE>) {

	if ($line =~ m/pipelineCL.py/) {
		chomp($line);
		$command = $line;
		#print STDOUT $line."\n";
		readline(FILE);
		}
	
	if ($line =~ m/Informatics\ Report\ Version\:/) {
		chomp($line);
		$pipelineversion = (split "\:", $line)[1];
		#print STDOUT $line."\n";		
		readline(FILE);
		}

	if ($line =~ m/Reading\ molecule\ stats\ from\ (.*)\:$/) {
		chomp($line);
		my $source = (split "/", $1)[-1];
		push @statst, $source;
		readline(FILE);
		# Reading molecule stats from /home/bionano/data/2312_ArabicaWS-2.5.1/2312-Arabica_FM/output/all.bnx:
		# read until empty line then store
		my @data = ();

		while ( $line = <FILE> ) {
			chomp($line);
			my ($title, $value);
			last if $line =~ /^$/;
			my @fields =  split(":", $line);
			($title = $fields[0]) =~ s/^\s+|\s+$//g;
			($value = $fields[1]) =~ s/^\s+|\s+$//g;
			push @data, [( $title, $value )];
			}
		@data = transpose(\@data);
		push @stats, [@data];
		}

	if ($line =~ m/Stage\ Summary\:\ Characterize(.*\ refineFinal1)$/) {
		# last two blocks
		chomp($line);
		my $source = $1;
		push @finalt, $source;
		# Stage Summary: CharacterizeDefault|Final refineFinal1
		# read until empty line then store
		my @data = ();

		while ( $line = <FILE> ) {
			chomp($line);
			my ($title, $value);
			last if $line =~ /^$/;
			my @fields =  split(":", $line);
			($title = $fields[0]) =~ s/^\s+|\s+$//g;
			($value = $fields[1]) =~ s/^\s+|\s+$//g;
			push @data, [( $title, $value )];
			}
		@data = transpose(\@data);
		push @final, [@data];
	} 
	
	if ($line =~ m/Stage\ Summary\:\ CharacterizeDefault\ (.*)$/) {
		# other stages
		chomp($line);
		my $source = $1;
		push @stagest, $source;
		# Stage Summary: CharacterizeDefault XXXXX
		# read until empty line then store
		my @data = ();

		while ( $line = <FILE> ) {
			chomp($line);
			my ($title, $value);
			last if $line =~ /^$/;
			my @fields =  split(":", $line);
			($title = $fields[0]) =~ s/^\s+|\s+$//g;
			($value = $fields[1]) =~ s/^\s+|\s+$//g;
			push @data, [( $title, $value )];
			}
		@data = transpose(\@data);
		push @stages, [@data];
		}

	if ($line =~ m/Molecules\ Aligned\ to\ (.*)\:/) {
		chomp($line);
		my $source = $1;
		push @alignst, $source;
		# Alignment summaries
		# read until empty line then store
		my @data = ();

		while ( $line = <FILE> ) {
			chomp($line);
			my ($title, $value);
			last if $line =~ /^$/;
			my @fields =  split(":", $line);
			($title = $fields[0]) =~ s/^\s+|\s+$//g;
			($value = $fields[1]) =~ s/^\s+|\s+$//g;
			push @data, [( $title, $value )];
			}
		@data = transpose(\@data);
		push @aligns, [@data];
		}

	if ($line =~ m/SV\ detect\:\ svdetect_exp_refineFinal1_sv/) {
		readline(FILE);
		my @data = ();
		while ( $line = <FILE> ) {
		chomp($line);
		last if $line =~ /^$/;
		push @data, $line;
		}
		push @svres, @data;
	}
}

#print Dumper \@statst;
#print Dumper \@stats;
#print Dumper \@stagest;
#print Dumper \@stages;
#print Dumper \@alignst;
#print Dumper \@aligns;
#print Dumper \@svres;

close FILE;

# print results to OUT
open OUT, "> $outfile" || die $!;

# print command and version
print OUT "## De Novo Assembly results\n";
print OUT "\n# Assembly Command :".$command."\n";
print OUT "\n# Pipeline version :".$pipelineversion."\n";


# print stats results
print OUT "\n# Molecule Alignment Results\n";
print OUT join("\t", "stats", @{@{$stats[0]}[0]}, "cvg (x)") . "\n";
for (my $idx=0; $idx < scalar @stats; $idx++) {
    print OUT join("\n", join("\t", $statst[$idx], @{@{$stats[$idx]}[1]})) . "\n";
}

# print assembly stages results
print OUT "\n# Assembly Stage Results\n";
print OUT join("\t", "assembly stages", @{@{$stages[0]}[0]}) . "\n";
for (my $idx=0; $idx < scalar @stages; $idx++) {
    print OUT join("\n", join("\t", $stagest[$idx], @{@{$stages[$idx]}[1]})) . "\n";
}

# print final stage results
print OUT "\n# Final Assembly Stage Results\n";
print OUT join("\t", "Final assembly stages", @{@{$final[0]}[0]}) . "\n";
for (my $idx=0; $idx < scalar @final; $idx++) {
    print OUT join("\n", join("\t", $finalt[$idx], @{@{$final[$idx]}[1]})) . "\n";
}

# print final stats results
print OUT "\n# Molecule alignments\n";
print OUT join("\t", "alignments", @{@{$aligns[0]}[0]}) . "\n";
for (my $idx=0; $idx < scalar @aligns; $idx++) {
    print OUT join("\n", join("\t", $alignst[$idx], @{@{$aligns[$idx]}[1]})) . "\n";
}

# print SV results
print OUT "\n# Structural differences with the provided Reference\n";
@svres = map {s/^\s+|\s+$//g; $_; } @svres;

for (my $idx=0; $idx < scalar @svres; $idx++) {
    print OUT join("\t", split(":", $svres[$idx])) . "\n";
}
close OUT;


####################################
########## example input ########### 

# /home/bionano/scripts//pipelineCL.py -U -d -T 228 -j 228 -N 4 -i 5 -a /home/bionano/data/2312_ArabicaWS-2.5.1/2312-Arabica_FM/optArguments_haplotype.xml -w -t /home/bionano/tools/ -l /home/bionano/data/2312_ArabicaWS-2.5.1/2312-Arabica_FM/output -b /home/bionano/data/2312_ArabicaWS-2.5.1/2312-Arabica_FM/Molecules.bnx -r /home/bionano/data/2312_ArabicaWS-2.5.1/2312-Arabica_FM/Ara-PB_v1_BspQI.cmap -C /home/bionano/scripts//clusterArguments.xml
# 
# Informatics Report Version: 2.3.0
# 
# Reading molecule stats from /home/bionano/data/2312_ArabicaWS-2.5.1/2312-Arabica_FM/output/all.bnx:
# Molecule Stats:
# N mols: 2266904
# Total len (Mb): 539117.200
# Avg len (kb)  :    237.821
# Mol N50 (kb)  :    234.416
# Lab (/100kb)  :      7.428
# 
# Stage Summary: CharacterizeDefault Assembly
# N Genome Maps: 5639
# Total Genome Map Len  (Mb): 2552.117
# Avg. Genome Map Len   (Mb):    0.453
# Median Genome Map Len (Mb):    0.303
# Genome Map n50        (Mb):    0.457
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    2.815
# N Genome Maps total align      : 2363 (0.42)
# Total Aligned Len (Mb)            : 1055.138
# Total Aligned Len / Ref Len       :    1.164
# Total Unique Aligned Len (Mb)     :  490.085
# Total Unique Aligned Len / Ref Len:    0.541
# 
# Reading molecule stats from /home/bionano/data/2312_ArabicaWS-2.5.1/2312-Arabica_FM/output/all_sorted.bnx:
# Molecule Stats:
# N mols: 2006989
# Total len (Mb): 487562.401
# Avg len (kb)  :    242.932
# Mol N50 (kb)  :    240.903
# Lab (/100kb)  :      7.840
# Ref    Cov (x):    537.794
# 
# Molecules Aligned to Reference:
# N mol align       :    554369
# Mol fraction align:         0.276
# Tot align len (Mb):     79379.6
# Effective Cov (x) :        87.203
# Avg align len (kb):       143.2
# Fraction align len:         0.163
# Tot confidence    :   7192201.6
# Avg confidence    :        13.0
# Avg FP(/100kb)    :         0.44
# Avg FP ratio      :         0.052
# Avg FN ratio      :         0.110
# Avg bpp           :       495.8
# Avg sf            :         0.293
# Avg sd            :        -0.041
# Avg sr            :         0.022
# 
# Stage Summary: CharacterizeDefault refineB1
# N Genome Maps: 5699
# Total Genome Map Len  (Mb): 2139.645
# Avg. Genome Map Len   (Mb):    0.375
# Median Genome Map Len (Mb):    0.236
# Genome Map n50        (Mb):    0.401
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    2.360
# N Genome Maps total align      : 3228 (0.57)
# Total Aligned Len (Mb)            : 1752.961
# Total Aligned Len / Ref Len       :    1.934
# Total Unique Aligned Len (Mb)     :  584.389
# Total Unique Aligned Len / Ref Len:    0.645
# 
# Stage Summary: CharacterizeDefault Merge0
# N Genome Maps: 2067
# Total Genome Map Len  (Mb): 1278.682
# Avg. Genome Map Len   (Mb):    0.619
# Median Genome Map Len (Mb):    0.306
# Genome Map n50        (Mb):    1.177
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.410
# N Genome Maps total align      : 1313 (0.64)
# Total Aligned Len (Mb)            : 1310.898
# Total Aligned Len / Ref Len       :    1.446
# Total Unique Aligned Len (Mb)     :  580.311
# Total Unique Aligned Len / Ref Len:    0.640
# 
# Stage Summary: CharacterizeDefault extension1_1
# N Genome Maps: 1748
# Total Genome Map Len  (Mb): 1290.375
# Avg. Genome Map Len   (Mb):    0.738
# Median Genome Map Len (Mb):    0.430
# Genome Map n50        (Mb):    1.184
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.423
# N Genome Maps total align      : 1205 (0.69)
# Total Aligned Len (Mb)            : 1429.205
# Total Aligned Len / Ref Len       :    1.576
# Total Unique Aligned Len (Mb)     :  583.633
# Total Unique Aligned Len / Ref Len:    0.644
# 
# Stage Summary: CharacterizeDefault Merge1
# N Genome Maps: 1413
# Total Genome Map Len  (Mb): 1205.436
# Avg. Genome Map Len   (Mb):    0.853
# Median Genome Map Len (Mb):    0.541
# Genome Map n50        (Mb):    1.324
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.330
# N Genome Maps total align      : 1066 (0.75)
# Total Aligned Len (Mb)            : 1271.101
# Total Aligned Len / Ref Len       :    1.402
# Total Unique Aligned Len (Mb)     :  583.915
# Total Unique Aligned Len / Ref Len:    0.644
# 
# Stage Summary: CharacterizeDefault extension1_2
# N Genome Maps: 1386
# Total Genome Map Len  (Mb): 1219.910
# Avg. Genome Map Len   (Mb):    0.880
# Median Genome Map Len (Mb):    0.568
# Genome Map n50        (Mb):    1.323
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.346
# N Genome Maps total align      : 1060 (0.76)
# Total Aligned Len (Mb)            : 1424.373
# Total Aligned Len / Ref Len       :    1.571
# Total Unique Aligned Len (Mb)     :  584.925
# Total Unique Aligned Len / Ref Len:    0.645
# 
# Stage Summary: CharacterizeDefault Merge2
# N Genome Maps: 1267
# Total Genome Map Len  (Mb): 1186.239
# Avg. Genome Map Len   (Mb):    0.936
# Median Genome Map Len (Mb):    0.616
# Genome Map n50        (Mb):    1.393
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.308
# N Genome Maps total align      : 1005 (0.79)
# Total Aligned Len (Mb)            : 1334.657
# Total Aligned Len / Ref Len       :    1.472
# Total Unique Aligned Len (Mb)     :  584.887
# Total Unique Aligned Len / Ref Len:    0.645
# 
# Stage Summary: CharacterizeDefault extension1_3
# N Genome Maps: 1252
# Total Genome Map Len  (Mb): 1191.209
# Avg. Genome Map Len   (Mb):    0.951
# Median Genome Map Len (Mb):    0.627
# Genome Map n50        (Mb):    1.393
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.314
# N Genome Maps total align      : 990 (0.79)
# Total Aligned Len (Mb)            : 1261.768
# Total Aligned Len / Ref Len       :    1.392
# Total Unique Aligned Len (Mb)     :  585.886
# Total Unique Aligned Len / Ref Len:    0.646
# 
# Stage Summary: CharacterizeDefault Merge3
# N Genome Maps: 1166
# Total Genome Map Len  (Mb): 1161.554
# Avg. Genome Map Len   (Mb):    0.996
# Median Genome Map Len (Mb):    0.692
# Genome Map n50        (Mb):    1.450
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.281
# N Genome Maps total align      : 938 (0.80)
# Total Aligned Len (Mb)            : 1018.932
# Total Aligned Len / Ref Len       :    1.124
# Total Unique Aligned Len (Mb)     :  585.235
# Total Unique Aligned Len / Ref Len:    0.646
# 
# Stage Summary: CharacterizeDefault extension1_4
# N Genome Maps: 1161
# Total Genome Map Len  (Mb): 1165.501
# Avg. Genome Map Len   (Mb):    1.004
# Median Genome Map Len (Mb):    0.697
# Genome Map n50        (Mb):    1.445
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.286
# N Genome Maps total align      : 939 (0.81)
# Total Aligned Len (Mb)            : 1083.045
# Total Aligned Len / Ref Len       :    1.195
# Total Unique Aligned Len (Mb)     :  585.734
# Total Unique Aligned Len / Ref Len:    0.646
# 
# Stage Summary: CharacterizeDefault Merge4
# N Genome Maps: 1135
# Total Genome Map Len  (Mb): 1156.211
# Avg. Genome Map Len   (Mb):    1.019
# Median Genome Map Len (Mb):    0.717
# Genome Map n50        (Mb):    1.463
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.275
# N Genome Maps total align      : 927 (0.82)
# Total Aligned Len (Mb)            : 1037.140
# Total Aligned Len / Ref Len       :    1.144
# Total Unique Aligned Len (Mb)     :  585.351
# Total Unique Aligned Len / Ref Len:    0.646
# 
# Stage Summary: CharacterizeDefault extension1_5
# N Genome Maps: 1132
# Total Genome Map Len  (Mb): 1158.927
# Avg. Genome Map Len   (Mb):    1.024
# Median Genome Map Len (Mb):    0.725
# Genome Map n50        (Mb):    1.459
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.278
# N Genome Maps total align      : 927 (0.82)
# Total Aligned Len (Mb)            : 1028.232
# Total Aligned Len / Ref Len       :    1.134
# Total Unique Aligned Len (Mb)     :  585.900
# Total Unique Aligned Len / Ref Len:    0.646
# 
# Stage Summary: CharacterizeDefault Merge5
# N Genome Maps: 1104
# Total Genome Map Len  (Mb): 1149.714
# Avg. Genome Map Len   (Mb):    1.041
# Median Genome Map Len (Mb):    0.745
# Genome Map n50        (Mb):    1.484
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    1.268
# N Genome Maps total align      : 908 (0.82)
# Total Aligned Len (Mb)            : 1002.889
# Total Aligned Len / Ref Len       :    1.106
# Total Unique Aligned Len (Mb)     :  585.144
# Total Unique Aligned Len / Ref Len:    0.645
# 
# Stage Summary: CharacterizeFinal refineFinal1
# Diploid N Genome Maps: 1933
# Diploid Genome Map Len        (Mb): 2163.044
# Diploid Avg. Genome Map Len   (Mb):    1.119
# Diploid Median Genome Map Len (Mb):    0.837
# Diploid Genome Map n50        (Mb):    1.557
# Haploid N Genome Maps: 1095
# Haploid Genome Map Len        (Mb): 1129.103
# Haploid Avg. Genome Map Len   (Mb):    1.031
# Haploid Median Genome Map Len (Mb):    0.731
# Haploid Genome Map n50        (Mb):    1.484
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    2.386
# N Genome Maps total align      : 1634 (0.85)
# Total Aligned Len (Mb)            :  681.562
# Total Aligned Len / Ref Len       :    0.752
# Total Unique Aligned Len (Mb)     :  301.471
# Total Unique Aligned Len / Ref Len:    0.333
# 
# Stage Summary: CharacterizeDefault refineFinal1
# Diploid N Genome Maps: 1933
# Diploid Genome Map Len        (Mb): 2163.044
# Diploid Avg. Genome Map Len   (Mb):    1.119
# Diploid Median Genome Map Len (Mb):    0.837
# Diploid Genome Map n50        (Mb):    1.557
# Haploid N Genome Maps: 1095
# Haploid Genome Map Len        (Mb): 1129.103
# Haploid Avg. Genome Map Len   (Mb):    1.031
# Haploid Median Genome Map Len (Mb):    0.731
# Haploid Genome Map n50        (Mb):    1.484
# Total Ref Len   (Mb):  906.597
# Total Genome Map Len / Ref Len :    2.386
# N Genome Maps total align      : 1633 (0.84)
# Total Aligned Len (Mb)            : 1701.781
# Total Aligned Len / Ref Len       :    1.877
# Total Unique Aligned Len (Mb)     :  588.191
# Total Unique Aligned Len / Ref Len:    0.649
# 
# SV detect: svdetect_exp_refineFinal1_sv
# map filename : len (Mb) : align len (Mb) (ratio) : N SV
#                  SV type :   N : N/tot
#                  complex : 9001 : 0.5839
#                 deletion :  261 : 0.0169
#            deletion_tiny :    2 : 0.0001
#                      end : 1522 : 0.0987
#                insertion :  679 : 0.0441
#           insertion_tiny :    1 : 0.0001
#                inversion :    2 : 0.0001
#        inversion_partial :    2 : 0.0001
#   trans_interchr_overlap :  116 : 0.0075
#    trans_interchr_repeat :  349 : 0.0226
#   translocation_interchr : 3479 : 0.2257
# N contigs  :      1
# N sv total :  15414  (15414.0000 per contig)
# N sv filter:   4423  (4423.0000 per contig)
# 
# Reading molecule stats from /home/bionano/data/2312_ArabicaWS-2.5.1/2312-Arabica_FM/output/all_sorted.bnx:
# Molecule Stats:
# N mols: 2006989
# Total len (Mb): 487562.401
# Avg len (kb)  :    242.932
# Mol N50 (kb)  :    240.903
# Lab (/100kb)  :      7.840
# Contig Cov (x):    225.406
# 
# Molecules Aligned to Assembly:
# N mol align       :    805554
# Mol fraction align:         0.401
# Tot align len (Mb):    164927.6
# Effective Cov (x) :        77.596
# Avg align len (kb):       204.7
# Fraction align len:         0.338
# Tot confidence    :  14078633.4
# Avg confidence    :        17.5
# Avg FP(/100kb)    :         0.62
# Avg FP ratio      :         0.062
# Avg FN ratio      :         0.112
# Avg bpp           :       501.9
# Avg sf            :         0.226
# Avg sd            :         0.004
# Avg sr            :         0.017
# 
