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
# Stephane Plaisance (VIB-NC+BITS) 2015/06/13; v2.00
# + add more counts
# + filter in steps and report all results
# Stephane Plaisance (VIB-NC+BITS) 2015/06/13; v2.01
# + fixed kb in size distribution
#
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Statistics::Descriptive;
use POSIX qw(strftime);

my $version = 2.01;
my $date = strftime "%m/%d/%Y", localtime;

# autoflush
$|=1;

############################
# handle command parameters
############################

getopts('i:l:x:s:p:m:h');
our($opt_i, $opt_l, $opt_x, $opt_s, $opt_p, $opt_m, $opt_h);

my $usage="You must provide a BNX file with -i
## Usage: bnxstats.pl <-i bnx-file>
# Additional optional parameters are:
# <-l minsize in kb (100)>
# <-x maxsize in kb (5000)>
# <-s minimal SNR (3.50)>
# <-p percentile (99)>
# <-m max-AvgIntensity (or percentile if undef)>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile = $opt_i || die $usage."\n";
my $minlen = $opt_l || 100;
my $maxlen = $opt_x || 5000;
my $minsnr = $opt_s || 3.5;
my $percentile = $opt_p || 99;
my $maxavgint = $opt_m || 0.6;
defined($opt_h) && die $usage."\n";

# open stream from BNX file
my $FILE = OpenArchiveFile ($inputfile) or die $!;
my $outpath = dirname($inputfile);

# remove possible suffixes from filename
my @sufx=( ".bnx", ".bnx.gzip", ".bnx.gz", ".bnx.zip" );
my $outbase = basename($inputfile, @sufx);

# include size limit, max intensity, and snr in file names
my $ai = defined($opt_m) ? $maxavgint : $percentile."pc-";
my $outfile = $outpath."/".$minlen."kb_".$ai."ai_".$minsnr."snr_".$outbase."_stats.txt";
open OUT, "> $outfile" || die $!;

# working variables
my $counttot = 0;
my $filttot = 0;
my $totnucl = 0;
my $filtnucl = 0;
my $first = 1;
my $stat = Statistics::Descriptive::Full->new();
my @BigArray = ();
our @result = ();
our $ailim;
our $spacer= join( "", "#", ("-" x 50 ), "\n");

################################
# parse data and store in array
################################

while (my $line = <$FILE>) {

	# check top line for "# BNX File Version:	1.2"
	if ($first == 1) {
		if ($line !~ /#\ BNX\ File\ Version:/) {
		die "$line\n This does not seem to be a bnx file";
			}
		$first = 0;
		}

	# header block
	if ($line =~ /^#/) {
		next;
		}

	# entering data part
	# load four lines in @data
	my @data = (); # store current four lines of data

	# test data consistency
	$line =~ /^0/ or die "aborted, does not seem to be a valid bnx format";
	$counttot++;
	push (@data, $line);

	## 0	1	94694.1 ...
	# 1	504.4	3008.2 ...
	# QX11	1.0645	1.3571 ...
	# QX12	0.0771	0.0778 ...

	# read three more lines
	for (my $d=1; $d<4; $d++) {
		$line = <$FILE> || die "premature end of file";
		push @data, $line;
		}

	# split and filter 0-line data
	my $zerol=$data[0];
	chomp($zerol);

	my (undef, undef, $Length, $AvgIntensity, $SNR, $NumberofLabels) =
		split(/\t/, $zerol);

	$totnucl += $Length;

	# where >= $minlen and >= $minsnr
	if ( $Length >= $minlen*1000 && $SNR >= $minsnr) {
		# add to filtered total length
		$filtnucl += $Length;
		$filttot++;
		}

	# store values in @LengthArray as (array of [arrays of 4 records])
	# [ length, labelfreq, SNR,avgInt ]
	my $labfreq = 100000*$NumberofLabels/$Length;

	# store results in @BigArray
	push (@BigArray,
		( [$Length, sprintf("%.2f", $labfreq), $SNR, $AvgIntensity] ) );
}
undef $FILE;

#################
# REPORT RESULTS
#################

# print title and header
my $topline = "# BioNanoGenomics molecule report (v".$version."), ".$date;
push (@result, $topline);
push (@result, "# input file: ".$outbase);

# print molecule counts in subsets
push (@result, "# tot-molecules: ".(format_num($counttot)));
push (@result, "# inputFile size (Gb): ".(sprintf("%.3f", $totnucl/1_000_000_000)) );
push (@result, "# Filters: min-length=".$minlen.
		"kb, max-len=".$maxlen.
		"kb, min SNR=".$minsnr.
		", max-AvgInt=".$maxavgint);

push (@result, $spacer);

my $bintitle = "\n# Merged Filtered Size Distribution (Size >= ".
	$minlen."kb, Size <= ".
	$maxlen."kb";

my @distheader=("molecules", "# of Molecules", "Quantity (Gb)",
	"Bin Mass Fraction(%)", "Labels(/100kb)"
	);

my @binheader = ("Length Bin", "# of Molecules", "Quantity (Gb)",
	"Bin Mass Fraction(%)", "Labels(/100kb)"
	);

my $quantheader=join("\t",
	(".", "min","25%","median","75%","max","mean", $percentile."%")
	);


############################
# Count in unfiltered data
############################

# reports stats for all molecules
report_stats("All molecules", \@BigArray);

###### subset filtered by size #####################

my @filtered = grep { $_->[0] >= $minlen*1000 &&
						$_->[0] <= $maxlen*1000 } @BigArray;

report_stats("length-filtered molecules", \@filtered);

###### subset filtered by size, SNR and avgInt######

my @fullfiltered = grep { $_->[2] >= $minsnr &&
						$_->[3] <= $ailim } @filtered;

report_stats("fully-filtered molecules", \@fullfiltered);

# output results
print OUT join("\n", @result)."\n";

close OUT;

exit 0;

##############
#### Subs ####

##############
sub OpenArchiveFile
{
    # $Filename passed in, handle to file passed out
    my $File = shift; # filename
    my $FH; # file handle

    if ($File =~ /.bnx$/) {
		open ($FH, "cat $File | ") or die ("$!: can't open file $File");
		} elsif ($File =~ /.bnx.zip$/) {
				open ($FH, "unzip -p $File | ") or die ("$!: can't open file $File");
			} elsif ($File =~ /(.bnx.gzip|.bnx.gz)$/) {
					open ($FH, "gzip -dc $File | ") or die ("$!: can't open file $File");
				} else {
						die ("$!: the file $File does seem to be a 'bnx' file");
					}
    return $FH;
}

##############
sub format_num {
	my $num = shift;
	$num = reverse $num;     # reverse the number's digits
	$num =~ s/(\d{3})/$1\'/g; # insert quote every 3 digits, from beginning
	$num = reverse $num;     # Reverse the result
	$num =~ s/^\'//;         # remove leading sep, if any
	return $num;
}

##############
sub bin_stats {
return undef unless ( scalar(@_) );

# count molecules
my $cnt = scalar(@_);

# count full bin size (Gb) [0]
my $binsum;
map { $binsum += $_->[0] } @_;
my $binsize = $binsum / 1_000_000_000;

# count mass fraction
my $massfr = $binsum / $totnucl;

# average label/100k [1]
my $sumlab;
map { $sumlab += $_->[1] } @_;
my $avglab = $sumlab/@_;

return ( $cnt,
	sprintf("%.3f", $binsize),
	sprintf("%.1f", 100*$massfr),
	sprintf("%.2f", $avglab)
	)
}

##############
sub get_stats {
# uses https://metacpan.org/module/Statistics::Descriptive
# test input
return undef unless ( scalar(@_) );

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@_);
my @result = ( $stat->quantile(0),
	$stat->quantile(1),
	$stat->quantile(2),
	$stat->quantile(3),
	$stat->quantile(4),
	$stat->mean(),
	scalar($stat->percentile($percentile))
	);

return @result;
}

##############
sub get_N50 {
return undef unless ( scalar(@_) );
my @sort = sort {$b <=> $a} @_;
my $totsum;
map { $totsum += $_ } @_;
my $cumsum; # cumulative sum
my $ranknum; # sorted rank
foreach my $curlength(@sort){
	$cumsum += $curlength;
	$ranknum++;
	if($cumsum >= $totsum/2){
		return ($curlength);
		last;
		}
	}
}

##############
sub report_stats {
	my($title, $array_ref) = @_;

	# global counts
	my $tot = scalar(@$array_ref);
	my $pc = sprintf("%.1f", 100*$tot/$counttot);

	# compute stats on columns 0 (Length) 3 (AvgIntensity), and 2 (AvgSNR)
	push (@result, "# distributions for ".$title." (N=".( format_num($tot) ).", ".$pc."%)");
	push (@result, $quantheader);

	# molecule lengths
	my @sizea =  map $_->[ 0 ],  @$array_ref;
	my $totnucl;
	map { $totnucl += $_ } @sizea;
	my @sizedist = map { sprintf("%d", $_/1000)} ( get_stats(@sizea) );
	push ( @result, join("\t", "length (kb) ", @sizedist) );

	# molecule averageIntensities
	my @aia = map $_->[ 3 ],  @$array_ref;
	my @aidist = map { sprintf("%.2f", $_)} ( get_stats(@aia) );
	push ( @result, join("\t", "AvgIntensity", @aidist) );

	# molecule SNR
	my @snra = map $_->[ 2 ],  @$array_ref;
	my @snrdist = map { sprintf("%.2f", $_)} ( get_stats(@snra) );
	push ( @result, join("\t", "AvgSNR      ", @snrdist) );

	# add lines to table for N50
	my $n50_length = get_N50(@sizea);
	push (@result, "\n# N50 molecule lengths (kb) : ".
		(sprintf("%.1f", $n50_length/1000)));

	# which of XXpc and maxai
	my $aipc = sprintf "%.1f", $aidist[6];
	$ailim = ( $aipc > $maxavgint ? $maxavgint : $aipc );
	push ( @result, "# refined max-AvgIntensity : ".$ailim );

	# Count by size-bins in unfiltered data
	my @bin20k = grep {$_->[0]>20000} @$array_ref;
	my @bin100k = grep {$_->[0]>100000} @$array_ref;
	my @bin150k = grep {$_->[0]>150000} @$array_ref;
	my @bin180k = grep {$_->[0]>180000} @$array_ref;
	my @bin250k = grep {$_->[0]>250000} @$array_ref;
	my @bin500k = grep {$_->[0]>500000} @$array_ref;
	my @bin150200k = grep {$_->[0]>=150000 && $_->[0]<199999} @$array_ref;

	# add lines to table for All molecules
	push (@result, "\n# size bin distribution");
	push (@result, join( "\t", @binheader) );
	push (@result, join( "\t", ( ">20k", bin_stats(@bin20k) ) ));
	push (@result, join( "\t", ( ">100k", bin_stats(@bin100k) ) ));
	push (@result, join( "\t", ( ">150k", bin_stats(@bin150k) ) ));
	push (@result, join( "\t", ( ">180k", bin_stats(@bin180k) ) ));
	push (@result, join( "\t", ( ">250k", bin_stats(@bin250k) ) ));
	push (@result, join( "\t", ( ">500k", bin_stats(@bin500k) ) ));
	push (@result, join( "\t", ( "150-200k", bin_stats(@bin150200k) ) ));

	# add spacer line
	push (@result, $spacer);

}
