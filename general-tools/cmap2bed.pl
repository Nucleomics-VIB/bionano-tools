#!/usr/bin/perl -w

# convert cmap to BED
# keep coverage data
#
## first version: 2014
## split BNX content to header file and four data files
## for use in [R]
## supports gzipped BNX files
# Stephane Plaisance (VIB-NC+BITS) 2015/12/02; v1.1
# added field-name as name in  BED
# visit our Git: https://github.com/BITS-VIB

# CMapId  ContigLength  NumSites  SiteID  LabelChannel  Position  StdDev  Coverage  Occurrence  GmeanSNR  lnSNRsd  SNR  count  SNR  ...
# 135     296576.6      30        1       1             20.0      440.1   7.0       7.0         12.9421   0.5523   0
# 135     296576.6      30        2       1             23509.0   250.4   8.0       5.0         10.0180   0.3294   0
# 135     296576.6      30        3       1             31570.9   362.1   10.0      10.0        14.5801   0.4625   0
# 135     296576.6      30        4       1             54253.6   298.5   10.0      7.0         12.7829   0.5828   0
# 135     296576.6      30        5       1             69351.2   235.6   11.0      10.0        12.2120   0.6260   0
# 135     296576.6      30        6       1             79377.4   154.1   13.0      12.0        20.1263   0.3291   0
# 135     296576.6      30        7       1             83950.3   288.0   13.0      11.0        11.7287   0.4744   0
# 135     296576.6      30        8       1             102420.2  152.8   18.0      17.0        10.0521   0.4836   0
# 135     296576.6      30        9       1             108167.1  207.0   20.0      18.0        29.3203   0.6851   0
# 135     296576.6      30        10      1             119554.1  312.2   21.0      19.0        20.6982   0.6213   0
# 135     296576.6      30        11      1             146719.3  176.0   21.0      17.0        12.1342   0.3440   0
# 135     296576.6      30        12      1             154667.2  110.5   21.0      15.0        10.2657   0.5366   0
# 135     296576.6      30        13      1             157192.9  152.9   21.0      11.0        6.7995    0.3945   0
# 135     296576.6      30        14      1             162947.3  205.3   21.0      18.0        11.6323   0.5859   0
# 135     296576.6      30        15      1             174126.2  152.8   21.0      18.0        11.5555   0.4048   0
# 135     296576.6      30        16      1             179872.5  320.7   21.0      19.0        14.6636   0.4581   0
# 135     296576.6      30        17      1             208594.1  130.1   21.0      20.0        10.4443   0.4377   0
# 135     296576.6      30        18      1             212488.5  159.3   21.0      20.0        17.9254   0.4843   0
# 135     296576.6      30        19      1             218820.9  158.7   20.0      17.0        7.8968    0.3425   0
# 135     296576.6      30        20      1             225102.1  191.5   19.0      17.0        9.7032    0.5308   0
# 135     296576.6      30        21      1             234699.2  212.2   17.0      8.0         7.5147    0.3437   0
# 135     296576.6      30        22      1             246712.4  220.8   17.0      16.0        11.7142   0.3384   0
# 135     296576.6      30        23      1             259797.4  222.7   17.0      16.0        13.3972   0.5904   0
# 135     296576.6      30        24      1             273125.6  91.4    16.0      15.0        23.9480   0.6105   0
# 135     296576.6      30        25      1             274542.4  164.6   15.0      15.0        16.7142   0.7249   0
# 135     296576.6      30        26      1             280886.0  126.9   14.0      14.0        11.5927   0.4991   0
# 135     296576.6      30        27      1             283960.1  119.2   13.0      12.0        15.8956   0.7651   0
# 135     296576.6      30        28      1             286296.9  145.6   8.0       8.0         20.2689   0.1993   0
# 135     296576.6      30        29      1             288358.3  220.3   8.0       8.0         10.2494   0.2746   0
# 135     296576.6      30        30      1             294373.8  0.0     8.0       8.0         15.5044   0.4464   0
# 135     296576.6      30        31      0             296576.6  0.0     1.0       1.0         0.0000    0.0000   0

# BED

# cmapIF start(position of previous line +1) end(position 6) filename + Coverage(8)
# cmap1 1 20 filename + 7.0

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

# handle command parameters
getopts('i:s:');
our($opt_i, $opt_s);

my $usage = "Aim: Convert cmap data to BED5. You must provide a cmap file with -i
# Usage: cmap2bed.pl <-i cmap-file> 
# Optional parameters:
# -s <field number for BED-score (0-based; default to coverage=7)>
# (v1.2) 1:ContigLength 2:NumSites 3:SiteID 4:LabelChannel 
#        5:Position\n#\t 6:StdDev 7:Coverage 8:Occurrence 
#        9:GmeanSNR 10:lnSNRsd";

my $inputfile = $opt_i || die $usage;
my $scorefield = $opt_s || 7;
grep( /^$scorefield$/, ( 1..10) ) || die "-s should be in [1..10]\n";

our %fieldnames = (
	1 => "ContigLength",
	2 => "NumSites",
	3 => "SiteID",
	4 => "LabelChannel",
	5 => "Position",
	6 => "StdDev",
	7 => "Coverage",
	8 => "Occurrence",
	9 => "GmeanSNR",
	10 => "lnSNRsd"
	);

# load full cmap data into an array(s)
open FILE, $inputfile or die $!;
my $outpath = dirname($inputfile);
my $outbase = basename($inputfile, ".cmap");
my $outfile = $outbase.".bed";
my $genome_file = $outbase."_genome.txt";

# result files
open OUT, "> $outfile" || die $!;

# declare variables
my $countcmap = 0;
our @comments = ();
our @header = ();
our @colnames = ();
our @coltypes = ();

# parse data file
while (my $line = <FILE>) {
	if ($line =~ /^#/) {
		parseheader($line);
		next;
		}
	# this is data
	$countcmap++;
 	# print first row
	my @field = split /\t/, $line;
	my $lblcnt = $field[2];
	my $start = 0;
	
	# print next rows until contig ends
	for (my $i=0; $i<$lblcnt; $i++) {
		print OUT join("\t ", $field[0], $start*1000, $field[5]*1000, $fieldnames{$scorefield}, $field[$scorefield])."\n";
		$start = $field[5];
		$line = <FILE>;
		@field = split /\t/, $line;
		}	
}
close FILE;
close OUT;

# print summary
print STDOUT "## input: ".basename($inputfile)."\n";
print STDOUT "## records: ".$countcmap."\n";
print STDOUT "## Colnames\n";
print STDOUT "# ".join(", ", @colnames)."\n";
# debug
print STDOUT "## Headers\n";
print STDOUT "# ".join("\n##", @header)."\n";
print STDOUT "## Comments\n";
print STDOUT "# ".join("\n##", @comments)."\n";

##############
#### Subs ####

sub parseheader {
my $line = shift;
chomp($line);
$line =~ s/\.\.\.$//;
# put header in array
my @arr = split("\t", $line);
if ($arr[0] =~ /^#h/) {
		@colnames = @arr[1 .. $#arr];
	} elsif ($arr[0] =~ /^#f/) {
		@coltypes = @arr[1 .. $#arr];
	} elsif ($#arr == 0) {
		push @comments, $line;
	} else {
		push @header, $line;
	}
}
