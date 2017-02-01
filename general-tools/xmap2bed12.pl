#!/usr/bin/perl -w

# xmap2bed12.pl (first version: 2016)
# convert xmap data to BED12 format
# report data in reference coordinates
# use aligned part of the reference for the thick line
# add unaligned lengths of the query cmap as thin lines to both end
#
# Stephane Plaisance (VIB-NC+BITS) 2016/05/26; v1.0
# Stephane Plaisance (VIB-NC+BITS) 2016/05/27; v1.1
# added translate BNG keys back to real names in Fasta assembly
# added custom colouring (https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7)
#
# visit our Git: https://github.com/Nucleomics-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use List::Util qw( min max );

# handle command parameters
getopts('i:x:k:r:vh');
our($opt_i, $opt_x, $opt_k, $opt_r, $opt_v, $opt_h);
our $version="1.1 (05-2016)";

my $usage = "Aim: Convert xmap data to BED12. You must provide a xmap file with -i
# script version:".$version."
# Usage: xmap2bed12.pl <-i xmap-file>
# Optional parameters (xmap v0.2) :
# -x <minimal value for score (default=0)>
# -r <RGB feature color 255,0,0=red (default=0 | black)>
# -k <key file (when provided, will rename the sequences to their original naming (default absent)>
# -v <report verbose summary>
# <-h to display this help>";

defined($opt_h) && die $usage . "\n";
my $inputfile = $opt_i || die $usage;
my $minscore = $opt_x || 0;
my $rgbcol = $opt_r || 0;
my $keyfile = $opt_k || undef;

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
	14 => "Alignment"
	);

###################################
# load full key data into an array

our %translate = ();

if (defined $keyfile) {
	open KEYS, $keyfile or die "# keyfile not found or not readable!";
	# store name translations to hash
	my $keycnt=0;
	defined($opt_v) && print STDERR "\n# loading key pairs\n";
	while (my $line = <KEYS>) {
		# ignore header lines and comments
		$line =~ s/\s+$//;
		next if ($line =~ /^#|^$|^CompntId/);
		# fill a hash with replacement numbers
		my @keys = split /\t/, $line;
		# CompntId	CompntName	CompntLength
		$translate{$keys[0]} = $keys[1];
		# print STDOUT $keys[0]." => ".$translate{$keys[0]}."\n";
		$keycnt++;
	}
	close KEYS;
	# test hash contains data
	$keycnt > 0 || die "# no data found in file!";
	defined($opt_v) && print STDERR "# loaded ".$keycnt." key rows\n";
}

# load xmap header and process content
open FILE, $inputfile or die $!;
my $outpath = dirname($inputfile);
my $outbase = basename($inputfile, ".xmap");
my $outfile = $outbase."_gt".$minscore."_xmap_bed12.bed";

# result files
open OUT, "> $outfile" || die $!;

# declare variables
my $countxmap = 0;
my $keptxmap = 0;
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
	# ignore empty lines
	next if ($line =~ /^\s*$/);
	# this is data
	$countxmap++;
	my @field = ( undef, (split /\t/, $line) );
	my $refid;
	if (defined ($keyfile)) {
		# translating from key
		$refid = $translate{$field[3]};
		} else {
			# keeping BNG key naming
			$refid = $field[3];
			}
	#my $recname=join("|", $field[1], $field[2], $field[10], $field[11], $field[12]);
	my $recname=join("|", $field[1], $field[2], $field[11], $field[12]);
	my $confid = $field[9];
	my $orient = $field[8];
	
	# print to BED 12 format
	# report the reference coordinates even when the aligned segments are unequal in size
	# name field= XmapEntryID|QryContigID|HitEnum|QryLen|RefLen|
	if ($field[9] > $minscore) {
		$keptxmap++;

		# forward cmaps
		if ($orient eq "+") {
			my $recstart=$field[6]-($field[4]-1);
			$recstart = ($recstart>0) ? $recstart : 0;  
			my $thickstart=$field[6];
			my $recend=$field[7]+($field[11]-$field[5]);
			$recend = ($recend>$field[12]) ? $field[12] : $recend; 
			my $thickend=$field[7];
			my $blocsize=$recend-$recstart;
			my $blockstart=0;
			print OUT join("\t", 
				$refid,
				int($recstart),
				int($recend),
				"\"".$recname."\"",
				$confid,
				$orient,
				int($thickstart),
				int($thickend),
				$rgbcol,
				"1",
				int($blocsize).",",
				int($blockstart).","
				)."\n";
			} else {
			# reversed cmaps (start and end coordinates are swapped ????)
			my $recstart=$field[6]-($field[11]-$field[4]);
			$recstart = ($recstart>0) ? $recstart : 0;  
			my $thickstart=$field[6];
			my $recend=$field[7]+($field[5]-1);
			$recend = ($recend>$field[12]) ? $field[12] : $recend; 
			my $thickend=$field[7];
			my $blocsizes=$recend-$recstart;
			my $blockstarts=0;
			print OUT join("\t",
				$refid,
				int($recstart),
				int($recend),
				"\"".$recname."\"",
				$confid,
				$orient,
				int($thickstart),
				int($thickend),
				$rgbcol,
				"1",
				int($blocsizes).",",
				int($blockstarts).","
				)."\n";			
		}
	}
}
close FILE;
close OUT;

# if -v verbose was requested
if (defined $opt_v) {
	# print summary
	print STDOUT "\n##### XMAP header information #####\n";
	print STDOUT "| input: ".basename($inputfile)."\n";
	print STDOUT "| records: ".$countxmap."\n";
	print STDOUT "| min-Confidence: ".$minscore."\n";
	print STDOUT "| kept-records: ".$keptxmap."\n";
	print STDOUT "| Colnames: ";
	print STDOUT join(", ", @colnames)."\n";
	# debug
	print STDOUT "\n##### Headers #####\n";
	print STDOUT "|".join("\n|", @header)."\n";
	print STDOUT "\n##### Comments #####\n";
	print STDOUT "|".join("\n|", @comments)."\n";
	print STDOUT "##############################\n";
	}

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
