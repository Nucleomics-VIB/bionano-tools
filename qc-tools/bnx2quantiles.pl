#!/usr/bin/perl -w

# bnx2quantiles.pl
# first version: 2016-05-31
# report data-distributions for several key BNX value types
# designed to work with BNX 1.2 format
#
# Stephane Plaisance (VIB-NC+BITS) 2016/06/01; v1.1
# added average inter-label distance
# Stephane Plaisance (VIB-NC+BITS) 2016/06/01; v1.2
# handle no label case leading to div/0 error (empty rows X11 and X12)
#
# visit our Git: https://github.com/Nucleomics-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Statistics::Descriptive;
use POSIX qw(strftime);

my $version = "1.2";
my $date = strftime "%m/%d/%Y", localtime;

# autoflush
$|=1;

############################
# handle command parameters
############################

getopts('i:p:P:h');
our($opt_i, $opt_p, $opt_P, $opt_h);

my $usage="You must provide a BNX file with -i
## Usage: bnx2quantiles.pl <-i bnx-file>
# script version:".$version."
# Additional optional parameters are:
# <-p additional low percentile (1)>
# <-P additional high percentile (99)>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile = $opt_i || die $usage."\n";
our $lowperc = $opt_p || 1;
our $highperc = $opt_P || 99;
defined($opt_h) && die $usage."\n";

# open stream from BNX file
my $FILE = OpenArchiveFile ($inputfile) or die $!;
my $outpath = dirname($inputfile);

# remove possible suffixes from filename
my @sufx=( ".bnx", ".bnx.gzip", ".bnx.gz", ".bnx.zip" );
my $outbase = basename($inputfile, @sufx);

my $outfile = $outpath."/".$outbase."_p-".$lowperc."_P-".$highperc."_quantiles.txt";
open OUT, "> $outfile" || die $!;

# working variables
my $counttot = 0;
my $totnucl = 0;
my $first = 1;
my @BigArray = ();
our @result = ();
our $spacer = join( "", "#", ("-" x 50 ), "\n");
our $cntln = 0;

my $quantheader = sprintf("%20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s",
	".", "min","25%","median","75%","max","mean", "stdev", "N50", $lowperc."%", $highperc."%", "skewness", "kurtosis"
	);

################################
# parse data and store in array
################################

while (my $line = <$FILE>) {
	
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

	######### analyze one record
	# 0-line data
	my $zerol=$data[0];
	chomp($zerol);

	my (undef, undef, $Length, $AvgIntensity, $SNR, $NumberofLabels) = split(/\t/, $zerol);
	$totnucl += $Length;
	my $labfreq = 100000*$NumberofLabels/$Length;
	
	# 1-line data
	my $onel=$data[1];
	chomp($onel);
	my (undef, @labpos) = split(/\t/, $onel);
	my @labdist = interlabel(\@labpos);
	my $avg_labdist = average(\@labdist);
	
	# X11 row data
	my $x11l = $data[2];
	chomp($x11l);
	my ( undef, @labsnr ) = split( /\t/, $x11l );
	# handle empty row (no labels on molecule)
	@labsnr || next;
	my $avg_labsnr = average(\@labsnr);
	
	# X12 row data
	my $x12l = $data[3];
	chomp($x12l);
	my ( undef, @labai ) = split( /\t/, $x12l );
	# handle empty row (no labels on molecule)
	@labai || next;
	my $avg_labai = average(\@labai);

	# store results in @BigArray
	push (@BigArray,
		( [$Length, 
		sprintf("%.2f",	$labfreq), 
		$SNR, 
		$AvgIntensity, 
		sprintf("%.4f",$avg_labsnr), 
		sprintf("%.4f",$avg_labai),
		sprintf("%.3f",$avg_labdist)] ) );
	
}
undef $FILE;

#################
# REPORT RESULTS
#################

# print title and header
my $topline = "# BioNanoGenomics data quantiles (v".$version."), ".$date;
push (@result, $topline);
push (@result, "# input file: ".basename($inputfile));

# print molecule counts in subsets
push (@result, "# tot-molecules: ".(format_num($counttot)));
push (@result, "# input file size (Gb): ".(sprintf("%.3f", $totnucl/1_000_000_000)) );
push (@result, "# additional percentiles:");
push (@result, "#   low_percentile: ".$lowperc);
push (@result, "#   high_percentile: ".$highperc);

push (@result, $spacer);

# reports stats for all molecules
report_stats("All labeled molecules", \@BigArray);

# output results to screen and file
print STDOUT join("\n", @result)."\n";
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
						die ("$!: the file $File does not seem to be a 'bnx' file");
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

###########
sub average {
	@_ == 1 or die ('Sub usage: $average = average(\@array);'); 
	my ($array_ref) = @_; 
	my $sum; 
	my $count = scalar @$array_ref; 

	foreach (@$array_ref) { $sum += $_; } 
	return $sum / $count; 
} 

###############
sub interlabel {
	@_ == 1 or die ('Sub usage: @interlabel = interlabel(\@array);');
	return undef unless ( scalar(@_) );
	my ($array_ref) = @_;
	my @result;
	# case only one label in current molecule
	(scalar @$array_ref) >2 || return 0;
	# shift first label out of array
	my $previous = shift @$array_ref;
	# remove last coordinate (= molecule end)
	splice @$array_ref, (scalar @$array_ref -1);
	# get each inter-label distance
	foreach (@$array_ref) {
		push @result, ($_-$previous);
		$previous=$_;
	}
	return(@result);
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
sub get_stats {
# uses https://metacpan.org/module/Statistics::Descriptive
# test input
return undef unless ( scalar(@_) );

my $n50 = get_N50(@_);
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@_);
my @result = ( 
	$stat->quantile(0),
	$stat->quantile(1),
	$stat->quantile(2),
	$stat->quantile(3),
	$stat->quantile(4),
	$stat->mean(),
	$stat->standard_deviation(),
	(defined( $n50 ) ? $n50 : "0"),
	(defined( scalar($stat->percentile($lowperc)) ) ? scalar($stat->percentile($lowperc)) : "0"),
	(defined( scalar($stat->percentile($highperc)) ) ? scalar($stat->percentile($highperc)) : "0"),
	(defined( $stat->skewness() ) ? $stat->skewness() : "0"),
	(defined( $stat->kurtosis() ) ? $stat->kurtosis() : "0")
	);

return @result;
}

##############
sub report_stats {
	my($title, $array_ref) = @_;

	# global counts
	my $tot = scalar(@$array_ref);
	my $pc = sprintf("%.1f", 100*$tot/$counttot);

	# compute stats on columns 0 (Length) 3 (AvgIntensity), 2 (AvgSNR), 
	push (@result, "# distributions for ".$title." (N=".( format_num($tot) ).", ".$pc."%)");
	push (@result, $quantheader);

	# molecule lengths [0]
	my @sizea =  map $_->[ 0 ],  @$array_ref;
	my $totnucl;
	map { $totnucl += $_ } @sizea;
	my @sizedist = map { sprintf("%d", $_/1000)} ( get_stats(@sizea) );
	push ( @result, sprintf("%20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s", 
		"length (kb)", @sizedist) );

	# molecule averageIntensities [3]
	my @aia = map $_->[ 3 ],  @$array_ref;
	my @aidist = map { sprintf("%.2f", $_)} ( get_stats(@aia) );
	push ( @result, sprintf("%20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s", 
		"molecAvgIntensity", @aidist) );

	# molecule SNR [2]
	my @snra = map $_->[ 2 ],  @$array_ref;
	my @snrdist = map { sprintf("%.2f", $_)} ( get_stats(@snra) );
	push ( @result, sprintf("%20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s", 
		"molecAvgSNR", @snrdist) );

	# molecule label frequency [1]
	my @labfr =  map $_->[ 1 ],  @$array_ref;
	my @labfrdist = map { sprintf("%.2f", $_)} ( get_stats(@labfr) );
	push ( @result, sprintf("%20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s", 
		"labelDensity", @labfrdist) );
	
	# label average distance [6]
	my @labdist =  map $_->[ 6 ],  @$array_ref;
	my @labdistdist = map { sprintf("%.2f", $_)} ( get_stats(@labdist) );
	push ( @result, sprintf("%20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s", 
		"labelAvgDistance", @labdistdist) );
		
	# label averageIntensities [5]
	my @laia = map $_->[ 5 ],  @$array_ref;
	my @laidist = map { sprintf("%.2f", $_)} ( get_stats(@laia) );
	push ( @result, sprintf("%20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s", 
		"labelAvgIntensity", @laidist) );

	# label snr [4]
	my @lsnra = map $_->[ 4 ],  @$array_ref;
	my @lsnrdist = map { sprintf("%.2f", $_)} ( get_stats(@lsnra) );
	push ( @result, sprintf("%20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s", 
		"labelAvgSNR", @lsnrdist) );
	
	# add spacer line
	push ( @result, $spacer );
	
	return @result;
}
