#!/usr/bin/perl -w

# mqr2bnx (first version: 06/2016)
# uses MoleculeQualityReport.xmap and associated RawMolecules.bnx data
# identifies the subset of molecules aligning to a reference
# saves the selected molecules to a new BNX file
#
# Stephane Plaisance (VIB-NC+BITS) 2016/06/12; v1.0
#
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

# autoflush
$|=1;

# handle command parameters
getopts('b:x:c:h');
our($opt_b, $opt_x, $opt_c, $opt_h);
our $version="1.0 (06-2016)";

my $usage = "Aim: Identifies molecules aligning to a reference from a MQR run and saves them to a new BNX file. You must provide the original BNX file with -b, the MQR xmap file with -x. Optionally, a minial confidence can be added to improve data quality.
# script version:".$version."
# Usage: mqr2bnx.pl <-b bnx-file> <-x xmap-file>
# Optional parameters:
# -c <minimal confidence score (default=0)>
# <-h to display this help>";

defined($opt_h) && die $usage . "\n";
my $bnxfile = $opt_b || die $usage;
my $xmapfile = $opt_x || die $usage;
my $minscore = $opt_c || "0";

# test input
grep( /^$minscore$/, ( 0..100) ) || die "-c should be in [0..100]\n";

# load xmap header and process content
open XMAP, $xmapfile or die $!;
open BNX, $bnxfile or die $!;

# name output
my $inpath = dirname($bnxfile);
my $inbase = basename($bnxfile, ".bnx");
my $outpath = $inpath."/".$inbase."_gt".$minscore."_xmap2bnx.bnx";

# result files
open OUT, "> $outpath" || die $!;

# xmap variables
my $countxmap = 0;
my %inxmap = ();

# parse xmap data file and collect BNX ID's
# 0=XmapEntryID 1=QryContigID 2=RefContigID 3=QryStartPos 4=QryEndPos
# 5=RefStartPos 6=RefEndPos 7=Orientation 8=Confidence 9=HitEnum

print STDERR "# parsing xmap and finding alignment records above confidence $minscore.\n";

while (my $line = <XMAP>) {
	# ignore header and empty lines
	next if ($line =~ /^#|^\s*$/);
	# this is data
	my @field = ( split /\t/, $line );
	if ($field[8] >= $minscore) {
		# will be extracted
		$countxmap++;
		# store BNX ID in hash
		$inxmap{$field[1]}++;
		}
	}
close XMAP;

# count aligning molecules
my @bnxids = keys(%inxmap);
my $uniquebnx = $#bnxids + 1;

print STDERR "# keeping $countxmap alignment records for $uniquebnx molecules.\n";

# parse BNX and keep header and selected molecules
my $count = 0;
my $kept = 0;
my $first = 1;

while ( my $line = <BNX> ) {
	# check top line for "# BNX File Version:	1.2"
	if ( $first == 1 ) {
		if ( $line !~ /#\ BNX\ File\ Version:/ ) {
			die "$line\n This does not seem to be a bnx file";
			}
		$first = 0;
		}

	# print header block
	if ( $line =~ /^#/ ) {
		print OUT $line;
		next;
		}

	# 0	1	94694.1 ...
	# 1	504.4	3008.2 ...
	# QX11	1.0645	1.3571 ...
	# QX12	0.0771	0.0778 ...
	#
	# 0=LabelChannel 1=MoleculeId 2=Length 3=AvgIntensity 4=SNR 5=NumberofLabels
	# 6=OriginalMoleculeId 7=ScanNumber 8=ScanDirection 9=ChipId 10=Flowcell

	# test data consistency
	$line =~ /^0/ or die "aborted, does not seem to be a valid bnx format";
	$count++;
	my @molecule = ();
	push @molecule, $line;

	# read three more lines
	for ( my $d = 1; $d < 4; $d++ ) {
		$line = <BNX> || die "premature end of file";
		push @molecule, $line;
	}

	# print to OUT when present in xmap
	my $bnxid = ( split( /\t/, $molecule[0]) )[1];
	if ( defined $inxmap{$bnxid} ) {
		$kept++;
  		map { print OUT "$_"; } @molecule;
  		}
	}

close BNX;
close OUT;

print STDERR "# saved $kept molecules from a total number of $count molecules in the input\n";
