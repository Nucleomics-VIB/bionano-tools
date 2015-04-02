#!/usr/bin/perl -w

# bnxfilter.pl
# first version: 2014-11-12
# filter a BioNanoGenomics RawMolecule file
# keep only molecules larger than minimum Length (100kb)
# keep only molecules with AvgIntensity smaller than max limit (0.4)
# keep only molecules with SNR greater than min limit (3.5)
# designed to work with BNX 1.2 format

# Stephane Plaisance (VIB-NC+BITS) 2015/04/02; v1.5
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

############################
# handle command parameters
############################

getopts('i:l:x:m:s:h');
our($opt_i, $opt_l, $opt_x, $opt_m, $opt_s, $opt_h);

my $usage="You must provide a BNX file with -i
## Usage: bnxfilter.pl <-i bnx-file>
# Additional optional parameters are:
# <-l min-length in kb (100)>
# <-x max-length in kb (5000)>
# <-m max-AvgIntensity (0.4)>
# <-s min-SNR (3.5)>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile = $opt_i || die $usage."\n";
my $minsize = $opt_l || 100;
my $minsizeb = $minsize * 1_000;
my $maxsize = $opt_x || 5000;
my $maxsizeb = $maxsize * 1_000;
my $maxavgint = $opt_m || 0.4;
my $minsnr = $opt_s || 3.5;
defined($opt_h) && die $usage."\n";

# open stream from BNX file
open FILE, $inputfile or die $!;
my $outpath = dirname($inputfile);
my $outbase = basename($inputfile);

# include size limit, max intensity, and snr in file names
my $outfile = $outpath."/filtered_".$minsize."kb_".
	$maxavgint."ai_".$minsnr."snr_".$outbase;
open OUT, "> $outfile" || die $!;

my $logfile = $inputfile."_".$minsize."kb_".
	${maxavgint}."ai_".$minsnr."snr_filtered-log.txt";
open OUT2, "> $logfile" || die $!;

# declare variables
my $counttot = 0;
my $goodcount = 0;
my $countshort = 0;
my $countlong = 0;
my $countlight = 0;
my $countlowsnr = 0;
my $first = 1;

################################
# parse data and store in array
################################

while (my $line = <FILE>) {

	# check top line for "# BNX File Version:	1.2"
	if ($first == 1) {
		if ($line !~ /#\ BNX\ File\ Version:/) {
		die "$line\n This does not seem to be a bnx file";
			}
		$first = 0;
		}

	# header block
	if ($line =~ /^#/) {
		# handle min-size change
		if ($line =~ /^# Min Molecule Length.*$/){
				print OUT "# Min Molecule Length (Kb):\t".$minsize."\n";
			} else {
				print OUT $line;
			}
		next;
		}

	# entering data part
	# load four lines in @data
	my @data=(); # store current four lines of data

	# test data consistency
	$line =~ /^0/ or die "aborted, does not seem to be a valid bnx format";
	$counttot++;
	push @data, $line;

	## 0	1	94694.1 ...
	# 1	504.4	3008.2 ...
	# QX11	1.0645	1.3571 ...
	# QX12	0.0771	0.0778 ...

	# read three more lines
	for (my $d=1; $d<4; $d++) {
		$line = <FILE> || die "premature end of file";
		push @data, $line;
		}

	# spit and filter 0-line data
	my $zerol=$data[0];
	chomp($zerol);
	my (undef, undef, $Length, $AvgIntensity, $snr, undef) =
		split(/\t/, $zerol);

	# filter
	if (($Length >= $minsizeb)
			&& ($Length <= $maxsizeb)
			&& ($AvgIntensity <= $maxavgint)
			&& ($snr >= $minsnr)) {
		# good record
		$goodcount++;
		foreach(@data) { print OUT "$_"; }
	} else {
		# too short
		if ($Length < $minsizeb) {
			$countshort++;
			}
		# too long
		if ($Length > $maxsizeb) {
			$countlong++;
			}
		# too fluorescent
		if ($AvgIntensity > $maxavgint) {
			$countlight++;
			}
		# low snr
		if ($snr < $minsnr) {
			$countlowsnr++;
			}
	}
}

close FILE;
close OUT;

##################
# prepare summary
##################

print OUT2 join("\t", "minLength", "AvgIntensity", "tot-molecules",
	"retained", "size-reject_short", "size-reject_long", "AvgInt-reject", "lowSNR-reject")."\n";

print OUT2 join("\t", $minsize, $maxavgint, $counttot, $goodcount,
	$countshort, $countlong, $countlight, $countlowsnr)."\n";

close OUT2;

exit 0;
