#!/usr/bin/perl -w

# labeldensity.pl
# Search for nicking enzyme sites in multifasta (reqired: restrict2bed.pl)
# create genome intervals from multifasta (reqired: fasta2chromsizes.pl)
# create windows (reqired: bedtools makewindows)
# compare both bed files and compute for each bin (reqired: bedtools map)
# report results in BED format visualisation
#
# Stephane Plaisance (VIB-NC+BITS) 2015/11/11; v1.00
#
# visit our Git: https://github.com/BITS-VIB

use warnings;
use strict;
use Getopt::Std;
use File::Basename;

# link the required scripts here
my $restrict2bed = "/work/MyScripts/restrict2bed.pl";
my $fasta2chromsizes = "/work/MyScripts/fasta2chromsizes.pl";
my $bedtools="/opt/biotools/bedtools/bin/bedtools";

############################
# handle command parameters
############################
getopts('i:n:l:b:h');
our ( $opt_i, $opt_n, $opt_l, $opt_b, $opt_h );

my $usage = "## Usage: labeldensity.pl <-i fasta-file>
# <-n 'nicker(s)', multiple allowed separated by ',')>
#  'Nt-BspQI' => 'GCTCTTC',
#  'Nt-BbvCI' => 'CCTCAGC',
#  'Nb-BsMI'  => 'GAATGC',
#  'Nb-BsrDI' => 'GCAATG'
# Additional optional parameters are:
# <-l minimal length for dna sequence (20000)>
# <-b bin width for computing label density (100000)>
# <-h to display this help>";

my $infile = $opt_i || die $usage . "\n";
my $nicker = $opt_n || die $usage . "\n";
my $minlen = $opt_l || 20000;
my $binwidth = $opt_b || 100000;
defined($opt_h) && die $usage . "\n";

# handle IO
my $inpath = dirname($infile);
my @sufx = ( ".fa", ".fasta", ".fsa" );
my $name = basename( $infile, @sufx );
my $outpath = $inpath."/".$name.".bed";
my $cmd;

# search nickers in fasta
$cmd="perl $restrict2bed -i $infile -l $minlen -n $nicker";
print STDERR "# ".(qq($cmd))."\n";
system($cmd);
print STDERR "\n\n";

my $nicking = $inpath."/".$name.".bed";

# create chromosome.length from multifasta
$cmd="perl $fasta2chromsizes -i $infile -l $minlen";
print STDERR "# ".(qq($cmd))."\n";
system($cmd);
print STDERR "\n\n";

my $chromsizes = $inpath."/".$name.".chrom.sizes";
my $windows = $inpath."/".$name."_".$binwidth."-bin.bed";

# create windows
$cmd="bedtools makewindows -g $chromsizes -w 100000 > $windows";
print STDERR "# ".(qq($cmd))."\n";
system($cmd);
print STDERR "\n\n";

my $result= $inpath."/".$name."_".$binwidth."-labeldensity.bed"

# compare and map
$cmd="perl $bedtools map -a $windows -b $nicking -c 5 -o sum > $result";
print STDERR "# ".(qq($cmd))."\n";
system($cmd);
print STDERR "\n\n";

exit 0;
