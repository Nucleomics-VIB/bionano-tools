#!/usr/bin/perl -w

## bnxsplitter.pl
## first version: 2015-12-01
## split BNX content to header file and four data files
## for use in [R]
## supports gzipped BNX files

# Stephane Plaisance (VIB-NC+BITS) 2015/12/01; v1.0
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

my $usage = "You must provide a BNX file with -i
## Usage: bnxsplitter.pl <-i bnx-file>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile = $opt_i || die $usage . "\n";
defined($opt_h) && die $usage . "\n";

# open stream from BNX file
# open FILE, $inputfile or die $!;
my $FILE = OpenArchiveFile($inputfile) or die $!;
my $outpath = dirname($inputfile);

# remove possible suffixes from filename
my @sufx = ( ".bnx", ".bnx.gzip", ".bnx.gz", ".bnx.zip" );
my $outbase = basename( $inputfile, @sufx );

# create output handles
open HEADER, "> $outpath/\.header.tsv" || die $!;
open ZERO, "> $outpath/\.zero.tsv" || die $!;
open ONE, "> $outpath/\.one.tsv" || die $!;
open QX11, "> $outpath/\.qx11.tsv" || die $!;
open QX12, "> $outpath/\.qx12.tsv" || die $!;

# declare variables
my $count = 0;
my $first = 1;

################################
# parse data and store in array
################################

while ( my $line = <$FILE> ) {

	# check top line for "# BNX File Version:	1.2"
	if ( $first == 1 ) {
		if ( $line !~ /#\ BNX\ File\ Version:/ ) {
			die "$line\n This does not seem to be a bnx file";
		}
		$first = 0;
	}

	# header block
	if ( $line =~ /^#/ ) {
		print HEADER $line;
		next;
	}
	
	## 0	1	94694.1 ...
	# 1	504.4	3008.2 ...
	# QX11	1.0645	1.3571 ...
	# QX12	0.0771	0.0778 ...

	# test data consistency
	$line =~ /^0/ or die "aborted, does not seem to be a valid bnx format";
	$count++;
	my @molecule = ();
	push @molecule, $line;

	# read three more lines
	for ( my $d = 1; $d < 4; $d++ ) {
		$line = <$FILE> || die "premature end of file";
		push @molecule, $line;
	}

	# split to four files
	print ZERO $molecule[0];
	print ONE $molecule[1];
	print QX11 $molecule[2];
	print QX12 $molecule[3];		
}	

# take care of handles neetly
undef $FILE;
close HEADER;
close ZERO;
close ONE;
close QX11;
close QX12;

exit 0;

##############
#### Subs ####

sub OpenArchiveFile {
	# $Filename passed in, handle to file passed out
	my $File = shift;	# filename
	my $FH;			  # file handle

	if ( $File =~ /.bnx$/ ) {
		open( $FH, "cat $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /.bnx.zip$/ ) {
		open( $FH, "unzip -p $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /(.bnx.gzip|.bnx.gz)$/ ) {
		open( $FH, "gzip -dc $File | " ) or die("$!: can't open file $File");
	} else {
		die("$!: the file $File does seem to be a 'bnx' file");
	}
	return $FH;
}
