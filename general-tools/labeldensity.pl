#!/usr/bin/perl -w

# labeldensity.pl
# Search for nicking enzyme sites in multifasta (reqired: restrict2bed.pl)
# create genome intervals from multifasta (reqired: fasta2chromsizes.pl)
# create windows (reqired: bedtools makewindows)
# compare both bed files and compute for each bin (reqired: bedtools map)
# sort BED files naturally requires a recent version of GNU sort
# report results in BED format visualisation
#
# Stephane Plaisance (VIB-NC+BITS) 2015/11/11; v1.00
# added direct support to create IGV track
# added return 0 when density is null (2016-04-19; v1.01)
#
# dependencies:
# restrict2bed.pl & fasta2chromsizes.pl (from: https://github.com/BITS-VIB/ngs-tools)
# bedtools (from: https://github.com/arq5x/bedtools2/releases)
# 2histo.R (https://github.com/BITS-VIB/plotting-tools)
#
# visit our Git: https://github.com/Nucleomics-VIB

use warnings;
use strict;
use Getopt::Std;
use File::Basename;

my $version="1.01";

# link the required scripts here
my $restrict2bed = `which restrict2bed.pl` || die "missing restrict2bed.pl, check your path\n";
chomp($restrict2bed);
my $fasta2chromsizes = `which fasta2chromsizes.pl` || die "missing fasta2chromsizes.pl, check your path\n";
chomp($fasta2chromsizes);
my $bedtools=`which bedtools` || die "missing bedtools, check your path\n";
chomp($bedtools);
my $plotR=`which 2histo.R` || die "missing 2histo.R, check your path\n";
chomp($plotR);

############################
# handle command parameters
############################
getopts('i:n:t:l:b:h');
our ( $opt_i, $opt_n, $opt_t, $opt_l, $opt_b, $opt_h );

my $usage = "## Usage: labeldensity.pl <-i fasta-file> <-n 'nicker(s)'>
# multiple allowed separated by ',')>
# eg. 'Nt-BspQI' => 'GCTCTTC',
# eg. 'Nt-BbvCI' => 'CCTCAGC',
# eg. 'Nb-BsMI'  => 'GAATGC',
# eg. 'Nb-BsrDI' => 'GCAATG',
# eg. 'Nb-BssSI' => 'CACGAG'
# Additional optional parameters are:
# <-t title ('label-density')>
# <-l minimal length for dna sequence (20000)>
# <-b bin width for computing label density (100000)>
# <-h to display this help>";

my $infile = $opt_i || die $usage . "\n";
my $nicker = $opt_n || die $usage . "\n";
my $title = $opt_t || "label-density";
my $minlen = $opt_l || 20000;
my $binwidth = $opt_b || 100000;
defined($opt_h) && die $usage . "\n";

# handle IO
my $inpath = dirname($infile);
my @sufx = ( ".fa", ".fasta", ".fsa" );
my $name = basename( $infile, @sufx );
my $cmd;

# search nickers in fasta
my $nicking = $inpath."/".$name."-".$title.".bed";

$cmd="perl $restrict2bed -i $infile -l $minlen -n $nicker | \
	sort -k 1V,1 -k 2n,2 -k 3n,3 > $nicking";
print STDERR "# ".(qq($cmd))."\n";
system($cmd) && die "! failed creating nicking BED file";
print STDERR "\n\n";

# create chromosome.length from multifasta
my $chromsizes = $inpath."/".$name.".chrom.sizes";

$cmd="perl $fasta2chromsizes -i $infile -l $minlen | \
	sort -k 1V,1 -k 2n,2 > $chromsizes";
print STDERR "# ".(qq($cmd))."\n";
system($cmd) && die "! failed creating chrom.sizes from fasta";
print STDERR "\n\n";

# create windows
my $windows = $inpath."/".$name."_".$binwidth."-bin.bed";

$cmd="$bedtools makewindows -g $chromsizes -w $binwidth | \
	sort -k 1V,1 -k 2n,2 -k 3n,3 > $windows";
print STDERR "# ".(qq($cmd))."\n";
system($cmd) && die "! failed creating windows from chrom.sizes";
print STDERR "\n\n";

# compare and map; create two outputs
my $result=$inpath."/".$name."_".$binwidth."-".$title."-labeldensity.bed";

$cmd="$bedtools map -nonamecheck -a $windows -b $nicking -c 5 -o sum -null 0 | sort -k 1V,1 -k 2n,2 -k 3n,3 > $result";
print STDERR "# ".(qq($cmd))."\n";
system($cmd) && die "! failed summarizing nicking data in window-bins";
print STDERR "\n\n";

# plot density distribution using r
my $plot=$inpath."/".$name."_".$binwidth."-".$title."-labeldensity.png";

$cmd="(cut -f 4 $result | 2histo.R) && mv 2histo.png $plot";
system($cmd) && die "! failed plotting density counts";
print "# density counts were plotted to $plot\n\n";

# report density counts
my $density = $inpath."/".$name."_".$binwidth."-".$title."-labeldensity_counts.txt";

$cmd="cut -f 4 $result | sort | uniq -c | \
	sed -e 's/ *//' -e 's/ /\t/'| \
	awk -v w=$binwidth 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"; 
	print \"# density in \"w\"b stepping windows\"}{print \$2, \$1}' | \
	sort -n > $density";
system($cmd) && die "! failed reporting density counts";
print "# density counts were stored in $density\n\n";

# create IGV track version
my $igv=$inpath."/".$name."_".$binwidth."-".$title."-labeldensity.igv";

$cmd="echo \"#track name=$title-$binwidth-density\" > $igv; awk -v title=$title 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"} {print \$1,\$2,\$3,title,\$4}' $result >> $igv";
print STDERR "# ".(qq($cmd))."\n";
system($cmd) && die "! failed creating IGV track";
print STDERR "\n\n";

exit 0;
