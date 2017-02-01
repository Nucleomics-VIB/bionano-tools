#!/usr/bin/perl -w

# cmap2bed.pl (first version: 2014)
# convert cmap to BED5 (strand is set to '.')
# keep coverage value or other <user-selected> field as BED 'score'
#
# Stephane Plaisance (VIB-NC+BITS) 2015/12/02; v1.3
# added field-name as name in BED
#
# visit our Git: https://github.com/Nucleomics-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

# handle command parameters
getopts('i:s:h');
our($opt_i, $opt_s, $opt_h);

my $usage = "Aim: Convert cmap data to BED5. You must provide a cmap file with -i
# Usage: cmap2bed.pl <-i cmap-file> 
# Optional parameters:
# -s <field number for BED-score (0-based; default to coverage=7)>
# (v1.2) 1:ContigLength 2:NumSites 3:SiteID 4:LabelChannel 
#        5:Position\n#\t 6:StdDev 7:Coverage 8:Occurrence 
#        9:GmeanSNR 10:lnSNRsd
# <-h to display this help>";

defined($opt_h) && die $usage . "\n";
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
		print OUT join("\t ", $field[0], 
		$start*1000, 
		$field[5]*1000, 
		$fieldnames{$scorefield}, 
		$field[$scorefield],
		".")."\n";
		$start = $field[5];
		$line = <FILE>;
		@field = split /\t/, $line;
		}	
}
close FILE;
close OUT;

# print summary
print STDOUT "\n##### CMAP header information\n";
print STDOUT "| input: ".basename($inputfile)."\n";
print STDOUT "| records: ".$countcmap."\n";
print STDOUT "| Colnames: ";
print STDOUT join(", ", @colnames)."\n";
# debug
print STDOUT "# Headers\n";
print STDOUT join("\n|", @header)."\n";
print STDOUT "# Comments\n";
print STDOUT join("\n|", @comments)."\n";
print STDOUT "##############################\n";

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
