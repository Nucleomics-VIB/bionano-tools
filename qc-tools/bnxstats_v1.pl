#!/usr/bin/perl -w

# bnxstats.pl
# first version: 2014-11-12
# added N50 2015-01-20
# added SNR quantiles
# replaced geometric mean by mean
# add params in file name
# parse a BNX file and compute simple stats
# designed to work with BNX 1.2 format

# Stephane Plaisance (VIB-NC+BITS) 2015/05/05; v1.51
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Statistics::Descriptive;
use POSIX qw(strftime);

my $version = 1.5;
my $date = strftime "%m/%d/%Y", localtime;

############################
# handle command parameters
############################

getopts('i:l:s:p:h');
our ( $opt_i, $opt_l, $opt_s, $opt_p, $opt_h );

my $usage = "You must provide a BNX file with -i
## Usage: bnxstats.pl <-i bnx-file>
# Additional optional parameters are:
# <-l minsize in kb (100)>
# <-s minimal SNR (3.50)>
# <-p percentile (99)>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile  = $opt_i || die $usage . "\n";
my $minlen     = $opt_l || 100;
my $minsnr     = $opt_s || 3.5;
my $percentile = $opt_p || 99;
defined($opt_h) && die $usage . "\n";

# open stream from BNX file
open FILE, $inputfile or die $!;
my $outpath = dirname($inputfile);
my $outbase = basename($inputfile);

# include size limit, max intensity, and snr in file names
my $outfile =
  $outpath . "/" . $minlen . "kb_" . $minsnr . "snr_" . $outbase . "_stats.txt";
open OUT, "> $outfile" || die $!;

# working variables
my $counttot = 0;
my $filttot  = 0;
my $totnucl  = 0;
my $filtnucl = 0;
my $first    = 1;
my $stat     = Statistics::Descriptive::Full->new();
my @BigArray = ();
my @result   = ();

################################
# parse data and store in array
################################

while ( my $line = <FILE> ) {

    # check top line for "# BNX File Version:	1.2"
    if ( $first == 1 ) {
        if ( $line !~ /#\ BNX\ File\ Version:/ ) {
            die "$line\n This does not seem to be a bnx file";
        }
        $first = 0;
    }

    # header block
    if ( $line =~ /^#/ ) {
        next;
    }

    # entering data part
    # load four lines in @data
    my @data = ();    # store current four lines of data

    # test data consistency
    $line =~ /^0/ or die "aborted, does not seem to be a valid bnx format";
    $counttot++;
    push( @data, $line );

    ## 0	1	94694.1 ...
    # 1	504.4	3008.2 ...
    # QX11	1.0645	1.3571 ...
    # QX12	0.0771	0.0778 ...

    # read three more lines
    for ( my $d = 1; $d < 4; $d++ ) {
        $line = <FILE> || die "premature end of file";
        push @data, $line;
    }

    # split and filter 0-line data
    my $zerol = $data[0];
    chomp($zerol);

    my ( undef, undef, $Length, $AvgIntensity, $SNR, $NumberofLabels ) =
      split( /\t/, $zerol );

    $totnucl += $Length;

    # where >= $minlen and >= $minsnr
    if ( $Length >= $minlen * 1000 && $SNR >= $minsnr ) {

        # add to filtered total length
        $filtnucl += $Length;
        $filttot++;
    }

    # store values in @LengthArray as (array of [arrays of 3 records])
    my $labfreq = 100000 * $NumberofLabels / $Length;
    push( @BigArray,
          ( [ $Length, sprintf( "%.2f", $labfreq ), $SNR, $AvgIntensity ] ) );
}
close FILE;

###################################
# Create filtered molecule subsets
###################################

my @filtered =
  grep { $_->[0] >= $minlen * 1000 && $_->[2] >= $minsnr } @BigArray;
my $filteredpc = sprintf( "%.1f", 100 * $filttot / $counttot );

# From bins of pre-defined length limits
my @bin20k  = grep { $_->[0] > 20000  && $_->[2] >= $minsnr } @BigArray;
my @bin100k = grep { $_->[0] > 100000 && $_->[2] >= $minsnr } @BigArray;
my @bin150k = grep { $_->[0] > 150000 && $_->[2] >= $minsnr } @BigArray;
my @bin180k = grep { $_->[0] > 180000 && $_->[2] >= $minsnr } @BigArray;
my @bin250k = grep { $_->[0] > 250000 && $_->[2] >= $minsnr } @BigArray;
my @bin500k = grep { $_->[0] > 500000 && $_->[2] >= $minsnr } @BigArray;
my @bin150200k =
  grep { $_->[0] >= 150000 && $_->[0] < 199999 && $_->[2] >= $minsnr }
  @BigArray;

##################
# prepare summary
##################

##########################
# print title and header
##########################

my $topline = "# BioNanoGenomics molecule report (v" . $version . "), " . $date;
push( @result, $topline );
push( @result, "# input file: " . $outbase );

# print molecule count
push( @result, "# tot-molecules: " . $counttot );
push( @result,
      "# inputFile size (Gb): "
        . ( sprintf( "%.3f", $totnucl / 1_000_000_000 ) ) );
push( @result, "# Filters: min-length=" . $minlen . "kb, min SNR=" . $minsnr );
push( @result,
      "# filtered-molecules: " . $filttot . " (" . $filteredpc . "%)" );
push( @result,
      "# filtered size (Gb): "
        . ( sprintf( "%.3f", $filtnucl / 1_000_000_000 ) ) );

my @distheader = (
                   "molecules",
                   "# of Molecules",
                   "Quantity (Gb)",
                   "Bin Mass Fraction(%)",
                   "Labels(/100kb)"
                 );

my $bintitle =
    "\n# Merged Filtered Size Distribution (Size >= "
  . $minlen
  . "kb, Label SNR >= "
  . $minsnr . ")";

push( @result, $bintitle );
push( @result, join( "\t", @distheader ) );

# add lines to table
push( @result, join( "\t", ( "Input File", bin_stats(@BigArray) ) ) );
push( @result, join( "\t", ( "Filtered",   bin_stats(@filtered) ) ) );

push( @result, "\n# SNR-filtered molecules in size bins" );
my @distheader2 = (
                    "Length Bin",
                    "# of Molecules",
                    "Quantity (Gb)",
                    "Bin Mass Fraction(%)",
                    "Labels(/100kb)"
                  );

push( @result, join( "\t", @distheader2 ) );

# add lines to table
push( @result, join( "\t", ( ">20k",     bin_stats(@bin20k) ) ) );
push( @result, join( "\t", ( ">100k",    bin_stats(@bin100k) ) ) );
push( @result, join( "\t", ( ">150k",    bin_stats(@bin150k) ) ) );
push( @result, join( "\t", ( ">180k",    bin_stats(@bin180k) ) ) );
push( @result, join( "\t", ( ">250k",    bin_stats(@bin250k) ) ) );
push( @result, join( "\t", ( ">500k",    bin_stats(@bin500k) ) ) );
push( @result, join( "\t", ( "150-200k", bin_stats(@bin150200k) ) ) );

################################
# Size-binned for ALL molecules
################################

# molecule lengths
my @sizea = map $_->[0], @BigArray;
my @sizedist = map { sprintf( "%d", $_ ) } ( get_stats(@sizea) );

# compute N50 for all molecules
my $n50_length_all = get_N50(@sizea);

# molecule averageIntensities
my @aia = map $_->[3], @BigArray;
my @aidist = map { sprintf( "%.2f", $_ ) } ( get_stats(@aia) );

# molecule averageSNR
my @snra = map $_->[2], @BigArray;
my @snrdist = map { sprintf( "%.2f", $_ ) } ( get_stats(@snra) );

push( @result, "\n# N50 value for ALL molecules (" . $counttot . ")" );
push( @result,
      "molecules-length-N50 (kb): "
        . ( sprintf( "%.1f", $n50_length_all / 1000 ) ) );

# compute stats on columns 0 (Length) 3 (AvgIntensity), and 2 (AvgSNR)
push( @result,
      "\n# quantile distribution for ALL molecules (" . $counttot . ")" );

my $header = "variable\t"
  . join( "\t",
          ( "min", "25%", "median", "75%", "max", "mean", $percentile . "%" ) );
push( @result, $header );

my $length = "length (kb) \t" . join( "\t", @sizedist );
push( @result, $length );
my $avgint = "AvgIntensity\t" . join( "\t", @aidist );
push( @result, $avgint );
my $avgsnr = "AvgSNR      \t" . join( "\t", @snrdist );
push( @result, $avgsnr );

#####################################
# Size-binned for filtered molecules
#####################################

# molecule lengths
my @sizeaf = map $_->[0], @filtered;
my @sizedistf = map { sprintf( "%d", $_ ) } ( get_stats(@sizeaf) );

# compute N50 for filtered molecules
my $n50_length_filt = get_N50(@sizeaf);

# molecule averageIntensities
my @aiaf = map $_->[3], @filtered;
my @aidistf = map { sprintf( "%.2f", $_ ) } ( get_stats(@aiaf) );

# molecule averageSNR
my @snraf = map $_->[2], @filtered;
my @snrdistf = map { sprintf( "%.2f", $_ ) } ( get_stats(@snraf) );

push( @result, "\n# N50 value for size-filtered molecules (" . $filttot . ")" );
push( @result,
      "molecules-length-N50 (kb): "
        . ( sprintf( "%.1f", $n50_length_filt / 1000 ) ) );

# compute stats on columns 0 (Length) 3 (AvgIntensity), and 2 (AvgSNR)
push( @result,
      "\n# quantile distribution for filtered molecules (" . $filttot . ")" );

push( @result, $header );

my $lengthf = "length (kb) \t" . join( "\t", @sizedistf );
push( @result, $lengthf );
my $avgintf = "AvgIntensity\t" . join( "\t", @aidistf );
push( @result, $avgintf );
my $avgsnrf = "AvgSNR      \t" . join( "\t", @snrdistf );
push( @result, $avgsnrf );

# print results to STDOUT and to file
print STDOUT join( "\n", @result ) . "\n";
print OUT join( "\n", @result ) . "\n";
close OUT;

exit 0;

##################
sub bin_stats {
    my @data = @_;

    # count molecules
    my $cnt = scalar(@data);

    # count full bin size (Gb) [0]
    my $binsum;
    map { $binsum += $_->[0] } @data;
    my $binsize = $binsum / 1_000_000_000;

    # count mass fraction
    my $massfr = $binsum / $totnucl;

    # average label/100k [1]
    my $sumlab;
    map { $sumlab += $_->[1] } @data;
    my $avglab = $sumlab / @data;

    return ( $cnt,
             sprintf( "%.3f", $binsize ),
             sprintf( "%.1f", 100 * $massfr ),
             sprintf( "%.2f", $avglab ) );
}

sub get_stats {

    # uses https://metacpan.org/module/Statistics::Descriptive
    # test input
    return undef unless ( scalar(@_) );

    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@_);
    my @result = (
                   $stat->quantile(0),
                   $stat->quantile(1),
                   $stat->quantile(2),
                   $stat->quantile(3),
                   $stat->quantile(4),
                   $stat->mean(),
                   scalar( $stat->percentile($percentile) )
                 );
    return @result;
}

sub get_N50 {
    my @sort = sort { $b <=> $a } @_;
    my $totsum;
    map { $totsum += $_ } @_;
    my $cumsum;     # cumulative sum
    my $ranknum;    # sorted rank
    foreach my $curlength (@sort) {
        $cumsum += $curlength;
        $ranknum++;
        if ( $cumsum >= $totsum / 2 ) {
            return ($curlength);
            last;
        }
    }
}
