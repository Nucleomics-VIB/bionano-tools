#!/usr/bin/perl -w

## xmap2ALLMAPSbed.pl
## convert HybridScaffold XMAP to BED6 for ALLMAPS
## requires the key to name conversion for RefContigID's

# Stephane Plaisance (VIB-NC+BITS) 2017/06/30; v1.0
# visit our Git: https://github.com/Nucleomics-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

############################
# handle command parameters
############################

getopts('x:k:h');
our ($opt_x, $opt_k, $opt_h);

my $usage = "Aim: Convert XMAP data to BED6 for use with ALLMAPS
## Usage: xmap2ALLMAPSbed.pl 
# <-x NGScontigs_HYBRID_SCAFFOLD.xmap_sorted.xmap> (required)
# <-k NGS-Reference-to-cmap_key.txt (required)
# <-h to display this help>";

####################
# declare variables
####################

my $xmapfile = $opt_x || die $usage . "\n";
my $keyfile = $opt_k || die $usage."\n";
defined($opt_h) && die $usage . "\n";

#########################
# load keys from keyfile 
#########################

# Used to translate IDs back to original Contig names
my @keys = ();
my %translate = ();

print STDOUT "\n# loading key pairs\n";
open KEYS, $keyfile or die $!;
my $keycnt = 0;
while (my $line = <KEYS>) {
	$line =~ s/\s+$//;
	next if ($line =~ /^#|^$|^CompntId/);
	$keycnt++;
	# fill a hash with replacement numbers
	my @keys = split /\t/, $line;
	$translate{$keys[0]} = $keys[1];
	print STDOUT $keys[0]." => ".$translate{$keys[0]}."\n";
}
print STDOUT "\n# $keycnt key pairs loaded\n";

close KEYS;

# create output handle
my $outpath = dirname($xmapfile);
# remove possible suffixes from filename
my @sufx = ( ".xmap" );
my $outbase = basename( $xmapfile, @sufx );

# open file handles
open XMAP, $xmapfile or die $!;
open BED, ">".$outpath."/".$outbase."_shortname.bed" || die $!;
open BED2, ">".$outpath."/".$outbase."_fullname.bed" || die $!;

# declare variables
my $count = 0;
my $first = 1;

## input XMAP ##
# XmapEntryID
# QryContigID => add 'Super-Scaffold_' => $optname
# RefContigID => $translate{RefContigID} => $ngsname
# QryStartPos => $optstart
# QryEndPos => $optend
# RefStartPos
# RefEndPos
# Orientation => $alistrand
# Confidence => $aliscore
# HitEnum
# QryLen
# RefLen
# LabelChannel
# Alignment

## OUTPUT BED6 ##
# col1 = $ngsname: the scaffold the optical map match maps to (or the edited scaffold name for BED2)
# col2 = $optstart: the start of the optical map match on the scaffold in col1
# col3 = $optend: the end of the optical map  match on the scaffold in col1 
# col4 = $optname: the name of the optical map match
# col5 = $aliscore: score
# col6 = $alistrand: strand

# parse XMAP file and export BED6
while ( my $line = <XMAP> ) {

	# pass through header block
	if ( $line =~ /^#/ ) {
		next;
	}
	
	chomp($line);
	
	# split XMAP row
	my ( 
		$XmapEntryID,
		$QryContigID,
		$RefContigID,
		$QryStartPos,
		$QryEndPos,
		$RefStartPos,
		$RefEndPos,
		$Orientation,
		$Confidence,
		$HitEnum,
		$QryLen,
		$RefLen,
		$LabelChannel,
		$Alignment 
		) = split("\t", $line);
	
	# check if translate $QryContigID to $ngsname works (wrong key file provided!)
	defined $translate{$RefContigID} || die "# Contig $RefContigID is not found in the provided key file!";
	
	# prepare BED with data as-is
	
	my ( 
		$ngsname, 
		$optstart,
		$optend,
		$optname,
		$aliscore,
		$alistrand 
		) = ( 
		$translate{$RefContigID},
		$RefStartPos,
		$RefEndPos,
		"Super-Scaffold_".$QryContigID,
		$Confidence,
		$Orientation
		);
	
	# write to bed
	print BED join("\t", 
		( 
		$ngsname, 
		$optstart,
		$optend,
		$optname,
		$aliscore,
		$alistrand
		) 
	)."\n";

	print BED2 join("\t", 
		( 
		$ngsname."_subseq_".$RefStartPos.":".$RefEndPos, 
		$optstart,
		$optend,
		$optname,
		$aliscore,
		$alistrand 
		)
	)."\n";
}

# close handles
close XMAP;
close BED;
close BED2;

exit 0;
