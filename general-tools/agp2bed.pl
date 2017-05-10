#!/usr/bin/perl -w

## agp2bed.pl
## first version: 2017-05-09
## extract original coordinate from AGP optical superscaffolds
## deduce true coordinates from '_subseq_start:end'

# Stephane Plaisance (VIB-NC+BITS) 2017/05/09; v1.0
# visit our Git: https://github.com/BITS-VIB

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

my $usage = "Aim: Convert AGP data to BED. You must provide a AGP file with -i
## Usage: agp2bed.pl <-i AGP-file>
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
open BED, ">".$outpath."/".$outbase.".bed" || die $!;

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
	
	# split line in elements
	chomp($line);
	my @field = split("\t", $line);
	
	# case N-gap
	if ($field[6] eq "scaffold"){
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
	}
	
	# write to bed
	print BED join("\t", ( (split(/;/, $seqname))[0], $start, $end, $info, $strand ))."\n";
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
