#!/usr/bin/perl -w

# xmapisec.pl (first version: 2016_06_17)
# create a molecule ID list from a BNX file
# parse two 'xmap' files obtained from the MQR alignment of this BNX data against two ref-cmaps
# identify molecules IDs present in either and both xmap files
# <optional> set a minimal confidence required to identify a match
# split the original BNX data in four new BNX subsets:
# 1) molecules aligning to ref-cmap#1
# 2) molecules aligning to ref-cmap#2
# 3) molecules aligning to both ref-cmaps
# 4) molecules that did not align to any ref-cmap
#
# Stephane Plaisance (VIB-NC+BITS) 2016/06/17; v1.0
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use List::Util qw( min max );

# disable buffering to get output during long process (loop)
$|=1;

# handle command parameters
getopts('i:a:b:n:c:zh');
our($opt_i, $opt_a, $opt_b, $opt_n, $opt_c, $opt_z, $opt_h);
our $version="1.0 (2016-06-17)";

my $usage = "Aim: Identify molecules specific to two ref-cmaps, ubiquitous, or not-aligning

# Usage: xmapisac.pl <-i bnx-file> <-a first-xmap-file> <-b 2nd-xmap-file>
# script version:".$version."
# Additional optional parameters are:
# -n <prefix for the output files> (default='isec_')>
# -c <minimal confidence score to be considered (default='undef')>
# -z zip the output to save space> (default OFF)>
# <-h to display this help>";

# 1 => "XmapEntryID",
# 2 => "QryContigID",
# 3 => "RefContigID",
# 4 => "QryStartPos",
# 5 => "QryEndPos",
# 6 => "RefStartPos",
# 7 => "RefEndPos",
# 8 => "Orientation",
# 9 => "Confidence",
# 10 => "HitEnum",
# 11 => "QryLen",
# 12 => "RefLen",
# 13 => "LabelChannel",
# 14 => "Alignment"

defined($opt_h) && die $usage . "\n";
my $inputfile = $opt_i || die $usage;
my $mapfile1 = $opt_a || die $usage;
my $mapfile2 = $opt_b || die $usage;
my $confcut = $opt_c || undef;
my $zipit = defined($opt_z) || undef;

# name output
my @suffixlist=(".bnx", ".bnx.gz", ".bnx.zip");
my($outbase, $outpath, $suffix) = fileparse($inputfile, @suffixlist);
if (defined ($opt_n)) {
	$outbase=$opt_n;
	}

# create 4 output file names
my $align1 = $outbase."_align-1.bnx";
my $align2 = $outbase."_align-2.bnx";
my $align3 = $outbase."_align-both.bnx";
my $align0 = $outbase."_align-none.bnx";

# also create 5 list files
my $bnxlst = $outbase."_bnx-idlist.txt";
my $lst1 = $outbase."_idlist-1.txt";
my $lst2 = $outbase."_idlist-2.txt";
my $lst3 = $outbase."_idlist-both.txt";
my $lst0 = $outbase."_idlist-none.txt";

# parse first xmap file
my $MAP1 = OpenMAP($mapfile1) or die $!;
print STDERR "# loading XMAP#1 aligned molecule IDs into hash#1\n";

# variables
my $first = 1;
my $cntln = 0;
my $cntali = 0;
my @fields = ();
my %in1hash = ();

while ( my $line = <$MAP1> ) {
	# count lines
	$cntln++;
	# check top lines for "# XMAP File Version:    0.2"
	if ($first == 1) {
		if ($line !~ /^#\ XMAP\ File\ Version:\t0.2$/) {
			if ($cntln > 20) {
				die "## This does not seem to be a xmap file";
				} else {
					next;
					}
			}
	$first = 0;
	next;
	}

	# ignore any other comment line
	next if ($line =~ /^#/);
	
    # test data consistency
    $line =~ /^[0-9]+.*/ or die "## does not seem to be a xmap file";
    $cntali++;
	@fields = split( /\t/, $line);

	# check confidence cutoff from $confcut
	if (! defined($confcut) || $fields[8] >= $confcut) {
		$in1hash{$fields[1]}++;
    	}
	}

print STDERR "## finished parsing XMAP#1 data for ".$cntali." alignments\n";
undef $MAP1;

# parse second xmap file
my $MAP2 = OpenMAP($mapfile2) or die $!;
print STDERR "# loading XMAP#2 aligned molecule IDs into hash#2\n";

# reuse $first $cntln $cntali @fields
$first = 1;
$cntln = 0;
$cntali = 0;
@fields = ();
my %in2hash = ();

while ( my $line = <$MAP2> ) {
	# count lines
	$cntln++;
	# check top lines for "# XMAP File Version:    0.2"
	if ($first == 1) {
		if ($line !~ /^#\ XMAP\ File\ Version:\t0.2$/) {
			if ($cntln > 100) {
				die "## This does not seem to be a xmap file";
				} else {
					next;
					}
			}
	$first = 0;
	next;
	}

	# ignore any other comment line
	next if ($line =~ /^#/);

    # test data consistency
    $line =~ /^[0-9]+.*/ or die "## does not seem to be a xmap file";
    $cntali++;
	@fields = split( /\t/, $line);

	# check confidence cutoff from $confcut
	if (! defined($confcut) || $fields[8] >= $confcut) {
		$in2hash{$fields[1]}++;
    	}
    }

print STDERR "## finished parsing XMAP#2 data for ".$cntali." alignments\n";
undef $MAP2;

# create file handles for saving data
if ( defined($zipit) ) {
	open ALI1, " | gzip -c > " . $align1 . ".gz" || die $!;
	open ALI2, " | gzip -c > " . $align2 . ".gz" || die $!;
	open ALI3, " | gzip -c > " . $align3 . ".gz" || die $!;
	open ALI0, " | gzip -c > " . $align0 . ".gz" || die $!;
	} else {
		open ALI1, "> $align1" || die $!;
		open ALI2, "> $align2" || die $!;
		open ALI3, "> $align3" || die $!;
		open ALI0, "> $align0" || die $!;
		}

# also create ID lists
open LST1, "> $lst1" || die $!;
open LST2, "> $lst2" || die $!;
open LST3, "> $lst3" || die $!;
open LST0, "> $lst0" || die $!;

# count matches
my ( $cnt1, $cnt2, $cntboth, $cntno ) = (0, 0, 0, 0);

# parse BNX data and split based on hashes
print STDERR "# parsing BNX data and saving to files  ";

my $BNX = OpenBNX($inputfile) or die $!;
open BNXLST, "> $bnxlst" || die $!;

# variables
$first = 1;
$cntln = 0;
$cntali = 0;
@fields = ();
my $countmol = 0;

while ( my $line = <$BNX> ) {
	# count lines
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

	# print header to all outputs
	if ( $line =~ /^#/ ) {
		print ALI1 $line;
		print ALI2 $line;
		print ALI3 $line;
		print ALI0 $line;
		next;
		}

	# test data consistency
    $line =~ /^0\b/ or die "## does not seem to be a valid bnx format";
    $countmol++;
    $countmol =~ /00000$/ && print STDERR ".";
	# this is the zero-line
	my $molid = (split(/\t/, $line))[1];
	print BNXLST $molid."\n";
	# add $line to @data
	my @data = ();
	push @data, $line;
	# read three more lines in to @data
	for ( my $d = 1; $d < 4; $d++ ) {
    	$line = <$BNX> || die "premature end of file";
        push @data, $line;
    	}

	# split according to hashes
	if ( ! exists($in1hash{$molid}) && ! exists($in2hash{$molid}) ) {
		# no match
		print ALI0 map { "$_" } @data;
		print LST0 $molid."\n";
		$cntno++;
		} elsif ( exists($in1hash{$molid}) && exists($in2hash{$molid}) ) {
			# matching both
			print ALI3 map { "$_" } @data;
			print LST3 $molid."\n";
			$cntboth++;
			} elsif ( exists($in1hash{$molid}) ) {
				# matching first
				print ALI1 map { "$_" } @data;
				print LST1 $molid."\n";
				$cnt1++;
				} elsif ( exists($in2hash{$molid}) ) {
					print ALI2 map { "$_" } @data;
					print LST2 $molid."\n";
					$cnt2++;
					} else {
						die "## What's up Dr: $molid?\n";
						}
	}

print STDERR "\n## finished saving $countmol molecules to BNX files\n\n";

close BNXLST;
close ALI1;
close ALI2;
close ALI3;
close ALI0;

# report results
print STDERR "# Results\tcounts\t%total\n";
print STDERR "# aligning to #map1 only:\t$cnt1\t".
	(sprintf("%d%%",100*$cnt1/$countmol))."\n";
print STDERR "# aligning to #map2 only:\t$cnt2\t".
	(sprintf("%d%%",100*$cnt2/$countmol))."\n";
print STDERR "# aligning to both maps:\t$cntboth\t".
	(sprintf("%d%%",100*$cntboth/$countmol))."\n";
print STDERR "# aligning to none:\t$cntno\t".
	(sprintf("%d%%",100*$cntno/$countmol))."\n";

##############
#### Subs ####

sub OpenBNX {
    # $Filename passed in, handle to file passed out
    my $File = shift;    # filename
    my $FH;              # file handle

    if ( $File =~ /\.bnx$/ ) {
        open( $FH, "cat $File | " ) or die("$!: can't open file $File");
    } elsif ( $File =~ /\.bnx\.zip$/ ) {
        open( $FH, "unzip -p $File | " ) or die("$!: can't open file $File");
    } elsif ( $File =~ /(\.bnx\.gzip|\.bnx\.gz)$/ ) {
        open( $FH, "gzip -dc $File | " ) or die("$!: can't open file $File");
    } else {
        die("$!: the file $File does not seem to be a 'bnx' file");
    }
    return $FH;
}

sub OpenMAP {
    # $Filename passed in, handle to file passed out
    my $File = shift;    # filename
    my $FH;              # file handle

    if ( $File =~ /\.xmap$/ ) {
        open( $FH, "cat $File | " ) or die("$!: can't open file $File");
    } elsif ( $File =~ /\.xmap\.zip$/ ) {
        open( $FH, "unzip -p $File | " ) or die("$!: can't open file $File");
    } elsif ( $File =~ /(\.xmap\.gzip|\.xmap\.gz)$/ ) {
        open( $FH, "gzip -dc $File | " ) or die("$!: can't open file $File");
    } else {
        die("$!: the file $File does not seem to be a 'xmap' file");
    }
    return $FH;
}

######### example MQR command
# RefAligner
# 	-f
# 	-ref .../reference.cmap
# 	-i .../Molecules.bnx
# 	-o .../MoleculeQualityReport
# 	-nosplit 2
# 	-BestRef 1
# 	-biaswt 0
# 	-Mfast 0
# 	-FP 1.5
# 	-sf 0.2
# 	-sd 0.0
# 	-A 5
# 	-outlier 1e-4
# 	-endoutlier 1e-3
# 	-S
# 	-1000
# 	-sr 0.04
# 	-resbias 5 64
# 	-maxmem 256
# 	-M 3
# 	-minlen 150
# 	-maxlen 2000
# 	-minsites 8
# 	-minSNR 3
# 	-MaxIntensity 0.6
# 	-T 1e-7
# 	-maxthreads 48
# 	-hashgen 5 3 2.4 1.5 0.05 5.0 1 1 2
# 	-hash
# 	-hashdelta 10
# 	-hashmaxmem 256
# 	-insertThreads 8
# 	-usecolor 1
# 	-stdout
# 	-stderr
