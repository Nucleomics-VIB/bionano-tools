#!/usr/bin/perl -w

## agp2ALLMAPSbed.pl
## convert HybridScaffold AGP to BED6 for ALLMAPS

# Stephane Plaisance (VIB-NC+BITS) 2017/06/30; v1.0
# visit our Git: https://github.com/Nucleomics-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

############################
# handle command parameters
############################

# disable buffering to get output during long process (loop)
$|=1; 

getopts('i:h');
our ( $opt_i, $opt_h );

my $usage = "Aim: Convert AGP data to ALLMAP BED. You must provide a AGP file with -i
## Usage: agp2ALLMAPSbed.pl <-i AGP-file>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile = $opt_i || die $usage . "\n";
defined($opt_h) && die $usage . "\n";

# open stream from AGP file
# open FILE, $inputfile or die $!;
my $FILE = OpenArchiveFile($inputfile) or die $!;
my $outpath = dirname($inputfile);

# remove possible suffixes from filename
my @sufx = ( ".agp", ".agp.gzip", ".agp.gz", ".agp.zip" );
my $outbase = basename( $inputfile, @sufx );

# create output handle
open BED, ">".$outpath."/".$outbase."_ALLMAPS.bed" || die $!;

# declare variables
my $count = 0;
my $first = 1;

################################
# parse data and store in array
################################

while ( my $line = <$FILE> ) {

	# check top line for "# BNX File Version:	1.2"
	if ( $first == 1 ) {
		if ( $line !~ /##agp-version/ ) {
			die "$line\n This does not seem to be a valid AGP file";
		}
		$first = 0;
	}

	# pass through header block
	if ( $line =~ /^#/ ) {
		next;
	}
	
	#############################
	# split AGP line in elements
	#############################
	# col0: Obj_Name / Super-Scaffold_2
	# col1: Obj_Start / 1
	# col2: Obj_End / 634629
	# col3: PartNum  / 1
	# col4: Compnt_Type / W
	# col5: CompntId_GapLength / 000000F_018_pilon
	# col6: CompntStart_GapType / 1
	# col7: CompntEnd_Linkage / 634629
	# col8: Orientation_LinkageEvidence / +
	
	chomp($line);
	my @field = split("\t", $line);
	my $optname = $field[0];

	# consider only if a Super_Scaffold entry
	$optname =~ m/Super-Scaffold_/ || next;

	# skip N-gaps
	if ($field[4] eq "N"){
		# this is a GAP
		next;
		}
		
	# parse fields
	my $CompntId = $field[5];
	# correct coordinates in case of subseq
	my $seqname = $CompntId;
	my $coordinates;
	my ( $start, $end, $strand ) = @field[6..8];
	if ($CompntId =~ m/(.*)(_subseq_)(.*)/g){
		$seqname = $1;
		$coordinates = $3;
		# update coordinates based on subseq info
		( $start, $end ) = ( split(/:/, $coordinates) );
	}
	my $score = 0;
	# handle dovetail names with ';" and keep only first name
	# print BED join("\t", ( (split(/;/, $seqname))[0], $start, $end, $optname, $score, $strand ))."\n";
	print BED join("\t", ( (split(/\ q/, $seqname))[0], $start, $end, $optname, $score, $strand ))."\n";
	# print BED join("\t", ( $seqname, $start, $end, $optname, $score, $strand ))."\n";
}	

# take care of handles neetly
undef $FILE;
close BED;

exit 0;

##############
#### Subs ####

sub OpenArchiveFile {
	# $Filename passed in, handle to file passed out
	my $File = shift;	# filename
	my $FH;			  # file handle

	if ( $File =~ /.agp$/ ) {
		open( $FH, "cat $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /.agp.zip$/ ) {
		open( $FH, "unzip -p $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /(.agp.gzip|.bnx.gz)$/ ) {
		open( $FH, "gzip -dc $File | " ) or die("$!: can't open file $File");
	} else {
		die("$!: the file $File does seem to be a 'agp' file");
	}
	return $FH;
}
