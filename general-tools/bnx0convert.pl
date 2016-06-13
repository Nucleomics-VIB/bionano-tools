#!/usr/bin/perl -w

## bnx0convert.pl
## first version: 2016-06-13
## convert original BNX v0.1 to current v1.2 format
# adds X11 and X12 rows and fill with arbitrary values

# Stephane Plaisance (VIB-NC+BITS) 2016/06/13; v1.0
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

getopts('i:o:h');
our ( $opt_i, $opt_o, $opt_h );
our $version="1.0 (2016-06-13)";

my $usage = "Aim: Reformat old BNX 0.1 format to current version 1.2. Arbitrary values are used for backbone SNR and average intensity. You must provide a BNX file with -i
# script version:".$version."
## Usage: bnx0convert.pl <-i bnx-file>
# <-o run name (default 'UnknownRun'>
# <-h to display this help>";

####################
# declare variables
####################

# default values for new header and data
our $nick1="";
my $n=defined($opt_o) ? $opt_o : "UnknownRun";
our $sfolder="\\".$n."\\Detect Molecules";
our $serial="Converted from v0.1 data";
our $rundate="1/1/2015 0:00:00";
our $filt="Static";
our $minsnr=0;
our $softv="";
our $pxscan=68819821;
our $stretch=0.85;
our $baseperpxl="500.00";
our $chipid="00000,00000,1/1/2015,000000000";
our $minlen=0;
our $runid=1;

# default values for updated data
our $defai="0.1000";
our $defSNR="10.0000";
our $flowcell=1;

my $inputfile = $opt_i || die $usage . "\n";
defined($opt_h) && die $usage . "\n";

# open stream from BNX file
# open FILE, $inputfile or die $!;
my $FILE = OpenArchiveFile($inputfile) or die $!;
my $outdir = dirname($inputfile);

# remove possible suffixes from filename
my @sufx = ( ".bnx", ".bnx.gzip", ".bnx.gz", ".bnx.zip" );
my $outbase = basename( $inputfile, @sufx );
my $outpath = $outdir."/".$outbase."_v1.2.bnx";

# create output handles
open OUT, "> $outpath" || die $!;

# declare variables
my $count = 0;
my $first = 1;

################################
# parse data and store in array
################################

while ( my $line = <$FILE> ) {

	# check top line for "# BNX File Version:	0.1"
	if ( $first == 1 ) {
		if ( $line !~ /#\ BNX\ File\ Version:\t0.1/ ) {
			die "$line\n This does not seem to be a bnx v0.1 file";
		}
		$first = 0;
		# make and print new header once
		print OUT join("\n", createHeader())."\n";
	}

	# avoid old header block
	next if ( $line =~ /^#/ );

	# test data consistency
	$line =~ /^0/ or die "aborted, does not seem to be a valid bnx v0.1 format";
	chomp($line);
	$count++;

	my @zerof = split("\t", $line);
	# check there is data there
	@zerof || die "no data found in #0f row!";

	# read the second line of data
	$line = <$FILE> || die "premature end of file";
	chomp($line);
	$line =~ /^1/ or die "aborted, does not seem to be a valid bnx v0.1 format";
	my @onef = split("\t", $line);
	@onef || die "no data found in #1f row!";

	# count labels for #0f
	my $labcnt = $#onef-1;

	# ignore molecule if no label
	if ( $labcnt<1 ) { print STDERR "# ".$labcnt." labels found \n"; next };

	####################
	# add missing data #
	####################

	# reformat #0f
	# Label Channel, MapID, Length
	# LabelChannel, MoleculeId, Length, AvgIntensity, SNR,
	#  NumberofLabels, OriginalMoleculeId, ScanNumber, ScanDirection, ChipId,
	#  Flowcell, RunId, GlobalScanNumber
	my @zeronew = ();
	push @zeronew, ( $zerof[0], $count,  $zerof[2], $defai, $defSNR, $labcnt,
	 $count, "1", "-1", $chipid, $flowcell, $runid, "1" );

	# add X11 and X12 lines and print out
	my @x11new = ( ($defSNR) x $labcnt );
	my @x12new = ( ($defai) x $labcnt );

	# print all out
	print OUT join("\t", @zeronew)."\n";
	print OUT join("\t", @onef)."\n";
	print OUT join("\t", "QX11", @x11new)."\n";
	print OUT join("\t", "QX12", @x12new)."\n";
	}

# take care of handles neetly
undef $FILE;
close OUT;

print STDERR "# converted $count BNX records from format 0.1 to format 1.2\n";

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

sub createHeader {
	# rebuild header for version 1.2
	my @header = ();

	push @header, "# BNX File Version:\t1.2";
	push @header, "# Label Channels:\t1";
	push @header, "# Nickase Recognition Site 1:\t".$nick1;
	push @header, "# Min Molecule Length (Kb):\t0";
	push @header, "# Label SNR Filter Type:\t".$filt;
	push @header, "# Min Label SNR:\t".$minsnr;
	push @header, "# Software Version:\t".$softv;
	push @header, ( join("\t", "#rh", "SourceFolder", "InstrumentSerial",
		"Time", "NanoChannelPixelsPerScan", "StretchFactor", "BasesPerPixel",
		"NumberofScans", "ChipId", "Flowcell", "LabelSNRFilterType",
		"MinMoleculeLength", "MinLabelSNR",
		"RunId") );
	push @header, ( join("\t", "# Run Data", $sfolder, $serial, $rundate,
		$pxscan, $stretch, $baseperpxl, "1", $chipid, "1",
		$filt, $minlen, $minsnr, $runid) );
	push @header, ( join("\t", "# Quality Score QX01:", "SNR") );
	push @header, ( join("\t", "# Quality Score QX02:", "Ave Intensity") );
	push @header, ( join("\t", "#0h", "LabelChannel", "MoleculeId", "Length",
		"AvgIntensity", "SNR", "NumberofLabels", "OriginalMoleculeId",
		"ScanNumber", "ScanDirection", "ChipId", "Flowcell", "RunId",
		"GlobalScanNumber") );
	push @header, ( join("\t", "#0f", "int", "int", "float", "float",
		"float", "int", "int", "int", "int",
		"string", "int", "int", "int") );
	push @header, ( join("\t", "#1h", "LabelChannel", "LabelPositions[N]") );
	push @header, ( join("\t", "#1f", "int", "float") );
	push @header, ( join("\t", "#Qh", "QualityScoreID", "QualityScores[N]") );
	push @header, ( join("\t", "#Qf", "str", "float") );
	return @header;
	}