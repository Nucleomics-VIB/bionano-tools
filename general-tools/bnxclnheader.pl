#!/usr/bin/perl -w

## bnxclnheader.pl
## first version: 2015-06-13
## clean header of a BioNanoGenomics RawMolecule file
## replace:
#   + replace [^a-zA-Z0-9.-] by '_'
## designed to work with BNX 1.2 format

# Stephane Plaisance (VIB-NC+BITS) 2015/06/13; v1.0
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

############################
# handle command parameters
############################

getopts('i:zh');
our ( $opt_i, $opt_z, $opt_h );

my $usage = "You must provide a BNX file with -i
## Usage: bnxclnheader.pl <-i bnx-file>
# Additional optional parameters are:
# <-z zip results (default OFF)>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile = $opt_i          || die $usage . "\n";
my $zipit     = defined($opt_z) || undef;
defined($opt_h) && die $usage . "\n";

# open stream from BNX file
# open FILE, $inputfile or die $!;
my $FILE = OpenArchiveFile($inputfile) or die $!;
my $outpath = dirname($inputfile);

# remove possible suffixes from filename
my @sufx = ( ".bnx", ".bnx.gzip", ".bnx.gz", ".bnx.zip" );
my $outbase = basename( $inputfile, @sufx );

# include size limit, max intensity, and snr in file names
my $outfile = $outpath . "/". $outbase . "_cleaned.bnx";

if ( defined($zipit) ) {
    open OUT, " | gzip -c > " . $outfile . ".gz" || die $!;
} else {
    open OUT, "> $outfile" || die $!;
}

my $logfile = $outpath . "/" . $outbase . "_cleaned-log.txt";

open OUT2, "> $logfile" || die $!;
print OUT2 "## cleaning the ".$inputfile." header\n#\n";

# declare variables
my $countrunlines = 0;
my $first         = 1;
my $smplname      = "";

################################
# parse data and store in array
################################

while ( my $line = <$FILE> ) {

    # check top line for "# BNX File Version:	1.2"
    if ( $first == 1 ) {
        if ( $line !~ /#\ BNX\ File\ Version:/ ) {
            die "$line\n This does not seem to be a bnx file";
        } else {
            $first = 0;
        }
    }

    # line to clean
    if ( $line =~ /^#\ Run\ Data/ ) {
        $countrunlines++;

        # convert to array
        my @fields       = split( /\t/, $line );
        my $sourceFolder = $fields[1];
        my @pathelement  = split( /\\/, $sourceFolder );
        my $count        = @pathelement;
        print OUT2 "# " . @pathelement[ ( $count - 2 ) ] . "\n";
        @pathelement[ ( $count - 2 ) ] =~ s/[^a-zA-Z0-9.-]/_/gi;
        print OUT2 "# -> " . @pathelement[ ( $count - 2 ) ] . "\n";
        $fields[1] = join("\\", @pathelement);
		print OUT join("\t", @fields);
		next;
    }

    # any other line
	print OUT $line;
}

print OUT2 "#\n# cleaned $countrunlines line(s)\n";
print STDOUT "# cleaned $countrunlines line(s)\n";

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
