#!/usr/bin/perl -w

## bnxreheader.pl
## first version: 2015-06-13
## replace header of a BioNanoGenomics RawMolecule file
## with valid header from a related file
## used to correct multiple syntax errors in files produced by IrysView filtering
## designed to work with BNX 1.2 format

# Stephane Plaisance (VIB-NC+BITS) 2016/02/26; v1.0
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

############################
# handle command parameters
############################

getopts('i:t:zh');
our ( $opt_i, $opt_t, $opt_z, $opt_h );

my $usage = "You must provide a BNX file with -i and a template with -t
## Usage: bnxreheader.pl <-i bnx-file> <-t template-file>
# Additional optional parameters are:
# <-z zip results (default OFF)>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile = $opt_i          || die $usage . "\n";
my $tmpltfile = $opt_t          || die $usage . "\n";
my $zipit     = defined($opt_z) || undef;
defined($opt_h) && die $usage . "\n";

# open stream from BNX file
# open FILE, $inputfile or die $!;
my $FILE = OpenArchiveFile($inputfile) or die $!;
my $outpath = dirname($inputfile);

# remove possible suffixes from filename
my @sufx = ( ".bnx", ".bnx.gzip", ".bnx.gz", ".bnx.zip" );
my $outbase = basename( $inputfile, @sufx );

# rename output
my $outfile = $outpath . "/". $outbase . "_header-swapped.bnx";

if ( defined($zipit) ) {
    open OUT, " | gzip -c > " . $outfile . ".gz" || die $!;
} else {
    open OUT, "> $outfile" || die $!;
}

my $logfile = $outpath . "/" . $outbase . "_cleaned-log.txt";
open OUT2, "> $logfile" || die $!;

#################################
# parse header from template file

print OUT2 "## reading header from ".$tmpltfile." \n#\n";

my $TMPL = OpenArchiveFile($tmpltfile) or die $!;
while ( my $line = <$TMPL> ) {
    # line to keep as header
    if ( $line !~ /^#/ ) {
        # stop when no # found
        last;
        } else {
        print OUT $line;
		next;
    }
}

#################################
# parse data from target file

print OUT2 "## reading data from ".$inputfile." \n#\n";

while ( my $line = <$FILE> ) {
	# ignore header lines
	next if $line =~ /^#/;

    # any other line
	print OUT $line;
}

close OUT;
close OUT2;

exit 0;

##############
#### Subs ####

sub OpenArchiveFile {
	# $Filename passed in, handle to file passed out
	my $File = shift;    # filename
	my $FH;              # file handle

	if ( $File =~ /.bnx$/ ) {
		open( $FH, "cat $File | " )
		  or die("$!: can't open file $File");
	} elsif ( $File =~ /.bnx.zip$/ ) {
		open( $FH, "unzip -p $File | " )
		  or die("$!: can't open file $File");
	} elsif ( $File =~ /(.bnx.gzip|.bnx.gz)$/ ) {
		open( $FH, "gzip -dc $File | " )
		  or die("$!: can't open file $File");
	} else {
		die("$!: the file $File does seem to be a 'bnx' file");
	}
	return $FH;
}
