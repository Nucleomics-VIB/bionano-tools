#!/usr/bin/perl -w

## dovetail2bed.pl
## first version: 2017-05-09
## extract original coordinate from Dovetail .table.txt & .input_breaks.txt files
## save in BED format

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

my $usage = "Aim: Convert Dovetail XXX.table.txt (and .input_breaks.txt from the same folder) to BED.
## Usage: dovetail2bed.pl <-i XXX.table.txt or XXX.input_breaks.txt>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile = $opt_i || die $usage . "\n";
defined($opt_h) && die $usage . "\n";

# open FILE, $inputfile or die $!;
my $outpath = dirname($inputfile);

# remove possible suffixes from filename
my @sufx = ( ".txt" );
my $outbase = basename( $inputfile, @sufx );
my $FILE = OpenArchiveFile($inputfile) or die $!;

# create output handle
open BED, ">".$outpath."/".$outbase.".bed" || die $!;

# declare variables
my $count = 0;
my $first = 1;

################################
# parse data and store in array
################################

while ( my $line = <$FILE> ) {

	# pass through header block
	if ( $line =~ /^#/ ) {
		next;
	}
	
	# split line in elements
	chomp($line);
	my @field = split("\t", $line);
	
	# parse fields
	my $info = join(";", @field);
	my $oriinfo = $field[1]."_subseq_".$field[2].":q".$field[3];
	my $seqname = $field[0];
	my ( $start, $end, $strand ) = @field[5,6,4];
	
	# write to bed
	print BED join("\t", ( (split(/;/, $seqname))[0], $start, $end, join("|", ($info, $oriinfo)), $strand ))."\n";
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

	if ( $File =~ /.txt$/ ) {
		open( $FH, "cat $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /.txt.zip$/ ) {
		open( $FH, "unzip -p $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /(.txt.gzip|.bnx.gz)$/ ) {
		open( $FH, "gzip -dc $File | " ) or die("$!: can't open file $File");
	} else {
		die("$!: the file $File does seem to be a 'txt' file");
	}
	return $FH;
}
