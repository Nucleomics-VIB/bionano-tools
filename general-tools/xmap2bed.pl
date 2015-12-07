#!/usr/bin/perl -w

# xmap2bed.pl (first version: 2014)
# convert the ref part of a xmap to BED5 (strand is set to '.')
# default use confidence value (9) for the score
# use aligned part of the reference for coordinates
#
# Stephane Plaisance (VIB-NC+BITS) 2015/12/04; v1.0
# added field-name as name in BED
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use List::Util qw( min max );

# handle command parameters
getopts('i:x:c:n:s:h');
our($opt_i, $opt_x, $opt_c, $opt_n, $opt_s, $opt_h);

my $usage = "Aim: Convert xmap data to BED5. You must provide a xmap file with -i
# Usage: xmap2bed.pl <-i xmap-file>
# Optional parameters (v0.2) :
# -x <minimal value for score (default=0)>
# -c <coordinate system used <'q'=query/'r'=ref> (default='r')
# -n <field number for BED-name (1-based; default to XmapEntryID=1)>
#        1:XmapEntryID 2:QryContigID 3:RefContigID 4:QryStartPos 5:QryEndPos
#        6:RefStartPos 7:RefEndPos 8:Orientation 9:Confidence
#       10: HitEnum 11:QryLen 12:RefLen 13:LabelChannel 14:Alignment
# -s <field number for BED-score (1-based; default to confidence=9)>
#        1:XmapEntryID 2:QryContigID 3:RefContigID 4:QryStartPos 5:QryEndPos
#        6:RefStartPos 7:RefEndPos 8:Orientation 9:Confidence
#       10: HitEnum 11:QryLen 12:RefLen 13:LabelChannel 14:Alignment
# <-h to display this help>";

defined($opt_h) && die $usage . "\n";
my $inputfile = $opt_i || die $usage;
my $minscore = $opt_x || 0;
my $coordinate = $opt_c || "r";
my $namefield = $opt_n || 1;
my $scorefield = $opt_s || 9;

# test input
grep( /^$coordinate$/, ( "q","r" ) ) || die "-c should be of 'q'/'r'\n";
grep( /^$namefield$/, ( 1..14) ) || die "-n should be in [1..14]\n";
grep( /^$scorefield$/, ( 1..14) ) || die "-s should be in [1..14]\n";

our %fieldnames = (
	1 => "XmapEntryID",
	2 => "QryContigID",
	3 => "RefContigID",
	4 => "QryStartPos",
	5 => "QryEndPos",
	6 => "RefStartPos",
	7 => "RefEndPos",
	8 => "Orientation",
	9 => "Confidence",
	10 => "HitEnum",
	11 => "QryLen",
	12 => "RefLen",
	13 => "LabelChannel",
	14 => "Alignment",
	);

# report choices
print STDOUT "\n##### BED-field options #####\n";
print STDOUT "| coordinates: ".$coordinate."\n";
print STDOUT "| seqlab: ".($coordinate eq 'q' ? $fieldnames{2} : $fieldnames{3})."\n";
print STDOUT "| start: ".($coordinate eq 'q' ? $fieldnames{4} : $fieldnames{6})."\n";
print STDOUT "| end: ".($coordinate eq 'q' ? $fieldnames{5} : $fieldnames{7})."\n";
print STDOUT "| name: ".$fieldnames{$namefield}."\n";
print STDOUT "| score: ".$fieldnames{$scorefield}."\n";
print STDOUT "| strand: ".$fieldnames{8}."\n";

# load xmap header and process content
open FILE, $inputfile or die $!;
my $outpath = dirname($inputfile);
my $outbase = basename($inputfile, ".xmap");
my $outfile = $outbase."_gt".$minscore.".bed";

# result files
open OUT, "> $outfile" || die $!;

# declare variables
my $countxmap = 0;
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
	$countxmap++;

	my @field = ( undef, (split /\t/, $line) );
	my $seqlab = ($coordinate eq 'q' ? $field[2] : $field[3]);
	my $start = ($coordinate eq 'q' ? $field[4] : $field[6]);
	my $end = ($coordinate eq 'q' ? $field[5] : $field[7]);
	my $coordstart = int( min($start, $end)+0.5);
	my $coordend = int( max($start, $end)+0.5);

	# print next rows until contig ends
	#
	# test $field[$scorefield]>$minscore and ?print
	if ($field[$scorefield] > $minscore) {
		print OUT join("\t", $seqlab,
			$coordstart,
			$coordend,
			$field[$namefield],
			$field[$scorefield],
			$field[8])."\n";
		}
}
close FILE;
close OUT;

# print summary
print STDOUT "\n##### XMAP header information #####\n";
print STDOUT "| input: ".basename($inputfile)."\n";
print STDOUT "| records: ".$countxmap."\n";
print STDOUT "| Colnames: ";
print STDOUT join(", ", @colnames)."\n";
# debug
print STDOUT "\n##### Headers #####\n";
print STDOUT "|".join("\n|", @header)."\n";
print STDOUT "\n##### Comments #####\n";
print STDOUT "|".join("\n|", @comments)."\n";
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
