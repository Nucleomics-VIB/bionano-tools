#!/usr/bin/perl -w

## agpFixCoords.pl
## first version: 2017-05-09
## Convert AGP coordinates when subseq encountered
## replace Bionano contig coordinates from '_subseq_start:end'

# Stephane Plaisance (VIB-NC+BITS) 2017/05/09; v1.0
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

my $usage = "Aim: Convert AGP coordinates when subseq encountered. You must provide a AGP file with -i
## Usage: agpFixCoords.pl <-i AGP-file>
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
open AGPOUT, ">".$outpath."/"."fixed_".$outbase.".agp" || die $!;

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
		print AGPOUT $line;
		next;
	}

	# pass through header block
	if ( $line =~ /^#/ ) {
		print AGPOUT $line;
		next;
	}
	
	#############################
	# split AGP line in elements
	#############################
	# col0: Obj_Name
	# col1: Obj_Start
	# col2: Obj_End 
	# col3: PartNum 
	# col4: Compnt_Type
	# col5: CompntId_GapLength
	# col6: CompntStart_GapType
	# col7: CompntEnd_Linkage
	# col8: Orientation_LinkageEvidence
	
	chomp($line);
	my @field = split("\t", $line);
	
	# skip N-gaps
	if ($field[5] eq "N"){
		# this is a GAP
		next;
		}
		
	# parse fields
	my $info = join(";", @field);
	my $CompntId = $field[5];
	my $seqname = $CompntId;
	my ( $start, $end, $strand ) = @field[6..8];
	my $coordinates;
	
	# correct coordinates in case of subseq
	if ($CompntId =~ m/(.*)(_subseq_)(.*)/g){
		$seqname = $1;
		$coordinates = $3;
		# update coordinates based on subseq info
		( $start, $end ) = ( split(/:/, $coordinates) );
		# update Contig name=nameOfField
		$seqname =~ s/_subseq_[[:digit:]]+\:[[:digit:]]+$//;
	}
	
	# write to AGP with correct coordinates
	print AGPOUT join("\t", 
		@field[0..4],
		$seqname,
		$start, 
		$end,
		$field[8]
		)."\n";
}	

# take care of handles neetly
undef $FILE;
close AGPOUT;

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
