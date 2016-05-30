#!/usr/bin/perl -w

## bnxfilter.pl
## first version: 2014-11-12
## filter a BioNanoGenomics RawMolecule file
## keep only molecules that:
#   + are larger than minimum Length (100kb)
#   + are smaller than maximum length (5000 kb)
#   + have at least N labels (6)
#   + have an Avg-Intensity smaller than max limit (0.6)
#   + have a SNR greater than min limit (3.5)
## report individual counts per type for failed molecules
## designed to work with BNX 1.2 format

# Stephane Plaisance (VIB-NC+BITS) 2015/04/02; v1.5
# Stephane Plaisance (VIB-NC+BITS) 2015/06/06; v1.6
# 	+ added support for archiving
# 	+ added filter min label count
# Stephane Plaisance (VIB-NC+BITS) 2016/05/25; v1.6.1
# + allow 10 additional header lines before '# BNX' 
#   when cli manipulations have added more comment rows

# Stephane Plaisance (VIB-NC+BITS) 2016/05/30; v2.0.0b
# + added filter label average-snr and label average-ai
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

############################
# handle command parameters
############################

getopts('i:l:x:m:s:n:zh');
our ( $opt_i, $opt_l, $opt_x, $opt_M, $opt_S, $opt_m, $opt_s, $opt_n, $opt_z, $opt_h );

my $usage = "You must provide a BNX file with -i
## Usage: bnxfilter.pl <-i bnx-file>
# Additional optional parameters are:
# <-l min-length in kb (100)>
# <-x max-length in kb (2500)>
# <-n min-nicks (6)>
# <-S min-molSNR (3.5)>
# <-M max-molAvgIntensity (0.6)>
# <-s min-labSNR (undef)>
# <-m max-labAvgIntensity (undef)>
# <-z zip results (default OFF)>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile = $opt_i          || die $usage . "\n";
my $minsize   = $opt_l          || 100;
my $minsizeb  = $minsize * 1_000;
my $maxsize   = $opt_x          || 2500;
my $maxsizeb  = $maxsize * 1_000;
my $minnck    = $opt_n          || 6;
my $minsnr    = $opt_S          || 3.5;
my $maxavgint = $opt_M          || 0.6;
my $minlsnr    = $opt_s         || undef;
my $maxlavgint = $opt_m         || undef;
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
my $outfile =
    $outpath
  . "/filtered_"
  . $minsize . "kb_"
  . $minnck . "nck_"
  . $minsnr . "SNR_"
  . $maxavgint . "AI_"
  . (defined($minlsnr) ? $minlsnr . "snr_" : "")
  . (defined($maxlavgint) ? $maxlavgint . "ai_" : "")
  . $outbase . ".bnx";

if ( defined($zipit) ) {
    open OUT, " | gzip -c > " . $outfile . ".gz" || die $!;
} else {
    open OUT, "> $outfile" || die $!;
}

my $logfile =
    $outpath . "/"
  . $outbase . "_"
  . $minsize . "kb_"
  . $minnck . "nck_"
  . $minsnr . "SNR_"
  . $maxavgint . "AI_"
  . (defined($minlsnr) ? $minlsnr . "snr_" : "")
  . (defined($maxlavgint) ? $maxlavgint . "ai_" : "")
  . "_filtered-log.txt";

open OUT2, "> $logfile" || die $!;

# declare variables
my $counttot    = 0;
my $goodcount   = 0;
my $countshort  = 0;
my $countlong   = 0;
my $countlight  = 0;
my $countlowsnr = 0;
my $countllight  = 0;
my $countllowsnr = 0;
my $countlownck = 0;
my $first       = 1;
our $cntln      = 0;

################################
# parse data and store in array
################################

while ( my $line = <$FILE> ) {

	# count header lines and abort if '# BNX' is not found in top 10 rows
	$cntln++;
	# check top line for "# BNX File Version:	1.2"
	if ($first == 1) {
		if ($line !~ /#\ BNX\ File\ Version:/) {
			if ($cntln > 10) {
				die "$line\n This does not seem to be a bnx file";
				} else {
					next;
				}
			}
		$first = 0;
		}

    # header block
    if ( $line =~ /^#/ ) {

        # handle min-size change
        if ( $line =~ /^# Min Molecule Length.*$/ ) {
            print OUT "# Min Molecule Length (Kb):\t" . $minsize . "\n";
        } else {
            print OUT $line;
        }
        next;
    }

    # entering data part
    # load four lines in @data
    my @data = ();    # store current four lines of data

    # test data consistency
    $line =~ /^0/ or die "aborted, does not seem to be a valid bnx format";
    $counttot++;
    push @data, $line;

    ## 0	1	94694.1 ...
    # 1	504.4	3008.2 ...
    # QX11	1.0645	1.3571 ...
    # QX12	0.0771	0.0778 ...

    # read three more lines
    for ( my $d = 1; $d < 4; $d++ ) {
        $line = <$FILE> || die "premature end of file";
        push @data, $line;
    }

    # spit and filter 0-line data
    my $zerol = $data[0];
    chomp($zerol);
    my ( undef, undef, $Length, $AvgIntensity, $snr, $nck, undef ) =
      split( /\t/, $zerol );

    # filter
    if (    ( $Length >= $minsizeb )
         && ( $Length <= $maxsizeb )
         && ( $AvgIntensity <= $maxavgint )
         && ( $snr >= $minsnr )
         && ( $nck >= $minnck ) )
    {
    
    	# test label values if cutoff are defined  
		# if cutoff for label snr, test average
		if (defined $minlsnr) {
		##############################
		# split and filter X11-line data
		my $x11l = $data[2];
		chomp($x11l);
		my ( undef, @labsnr ) = split( /\t/, $x11l );
		# take average
		my $avg_labsnr = map $_/@labsnr, @labsnr;
		# fail if less than cutoff
		if ($avg_labsnr < $minlsnr ) {
			$countllowsnr++;
			next
			}
		}
	
		# if cutoff for label averageIntensity, test average
		if (defined $maxlavgint) {
		##############################
		# split and filter X12-line data
		my $x12l = $data[3];
		chomp($x12l);
		my ( undef, @labai ) = split( /\t/, $x12l );
		my $avg_labai = map $_/@labai, @labai;
		# fail if more than cutoff
		if ($avg_labai > $maxlavgint) {
			$countllight++;
			next
			}
		}
	
        # good record
        $goodcount++;
        foreach (@data) { print OUT "$_"; }
    } else {

        # too short
        if ( $Length < $minsizeb ) {
            $countshort++;
        }

        # too long
        if ( $Length > $maxsizeb ) {
            $countlong++;
        }

        # too fluorescent
        if ( $AvgIntensity > $maxavgint ) {
            $countlight++;
        }

        # low snr
        if ( $snr < $minsnr ) {
            $countlowsnr++;
        }

        # low label frequency
        if ( $nck < $minnck ) {
            $countlownck++;
        }
    }
}

# take care of handles neetly
undef $FILE;
close OUT;

##################
# print summary
##################

print OUT2 join( "\t",
                 "minLength", "maxLength", "minLabelCount",
 "maxAvgIntensity", "minSNR", "tot-molecules",
 "retained", "size-reject_short", "size-reject_long",
 "lowNCK-reject", "mol_highAvgInt-reject", "mol_lowSNR-reject",
 "lab_lowAvgSNR-reject", "lab_highAvgInt-reject" )
  . "\n";
print OUT2 join( "\t",
                 $minsize, $maxsize, $minnck, $maxavgint,
  $minsnr, $counttot, $goodcount, $countshort,
  $countlong, $countlownck, $countlight, $countlowsnr,
  $countllowsnr, $countllight )
  . "\n";

close OUT2;

exit 0;

##############
#### Subs ####

sub OpenArchiveFile {

    # $Filename passed in, handle to file passed out
    my $File = shift;    # filename
    my $FH;              # file handle

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
