#!/usr/bin/perl

# Add N-regions's to both ends of each hybridScaffold fasta record
# based on the HYBRID_SCAFFOLD_trimHeadTailGap.coord file content
#
# Stephane Plaisance (VIB-NC+BITS) 2017/06/28; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

use warnings;
use strict;
use Getopt::Std;
use Bio::SeqIO;
use POSIX;

my $version = "1.0";

my $usage="## Usage: fastaPadHeadTrailGap.pl
## script version:".$version."
# <-i NGScontigs_HYBRID_SCAFFOLD.fasta file (required)>
# <-x NGScontigs_HYBRID_SCAFFOLD.xmap_sorted.xmap (required)>
# <-t NGScontigs_HYBRID_SCAFFOLD_trimHeadTailGap.coord file (required)>
# <-h to display this help>";

####################
# declare variables
####################
getopts('i:x:t:h');
our ($opt_i, $opt_x, $opt_t, $opt_h);
my $infile = $opt_i || die $usage."\n";
my $xmapfile = $opt_x || die $usage."\n";
my $gapfile = $opt_t || die $usage."\n";
defined($opt_h) && die $usage."\n";

#########################
# load HeadTrailGap hash 
#########################

#################
# (example head)
#################
# ##agp-version	2.0
# # Organism:   
# # Platform:     
# # Model:        
# # Enzyme(s):    
# # BioSample:    
# # BioProject:   
# Obj_Id	HeadTrimmedLength	TailTrimmedLength
# Super-Scaffold_15	6204	8855

my %headtrim =();
my %tailtrim = ();

print STDOUT "\n# loading key pairs\n";
open TRIMS, $gapfile or die $!;
my $trimcnt = 0;
while (my $line = <TRIMS>) {
	$line =~ s/\s+$//;
	next if ($line =~ /^#|^$|^Obj_Id/);
	$trimcnt++;
	# fill a hash with replacement numbers
	my @recs = split /\t/, $line;
	$headtrim{$recs[0]} = $recs[1];
	$tailtrim{$recs[0]} = $recs[2];
	printf STDOUT "%s\t%s\t%s\n", $recs[0], $recs[1], $recs[2];
}
print STDOUT "# => $trimcnt HeadTrimmedLength and TailTrimmedLength pairs loaded\n";
close TRIMS;

########################
# load xmap_sorted.xmap 
########################

# we need the xmap info to find reverse oriented Superscafolds for which head and trip gaps should be swapped
#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	\
#	Confidence	HitEnum	QryLen	RefLen	LabelChannel	Alignment

my %Orientation = ();
my %QryStartPos = ();

print STDOUT "\n# loading xmap info\n";
open XMAP, $xmapfile or die $!;
my $xmapcnt = 0;
while (my $line = <XMAP>) {
	$line =~ s/\s+$//;
	next if ($line =~ /^#|^$/);
	# fill a hashes with QryStartPos & Orientation
	my @recs = split /\t/, $line;
	if ( $recs[1] eq $recs[2] ) {
		$xmapcnt++;
		my $SID = "Super-Scaffold_". $recs[2];
		$QryStartPos{$SID} = floor($recs[4]);
		$Orientation{$SID} = $recs[7];
		printf STDOUT "%s:%s\t%s\t%s\n", $recs[2],  $SID, floor($recs[4]), $recs[7];
	}
}
print STDOUT "# => $xmapcnt SuperScaffold strands loaded\n";
close XMAP;

###########################################
# process NGScontigs_HYBRID_SCAFFOLD.fasta 
###########################################

my $seq_in = Bio::SeqIO -> new(-file => "$infile", -format => 'Fasta');
my $seq_out = Bio::SeqIO -> new(-file => "> untrimmed_$infile", -format => 'Fasta');

while ( my $seq = $seq_in->next_seq() ) {
	my $seqid = $seq->display_id();
	my $headt = 0;
	my $tailt = 0;
	print STDOUT "# processing \'$seqid\' : ";
	if ($Orientation{$seqid} eq "-") {
 		$tailt = $headtrim{$seqid};
 		#$headt = $tailtrim{$seqid} + $QryStartPos{$seqid};
 		$headt = $tailtrim{$seqid};
	} else {
		$headt = $headtrim{$seqid};
 		#$headt = $headtrim{$seqid} + $QryStartPos{$seqid};
 		$tailt = $tailtrim{$seqid};
	}
	printf STDOUT "%s:%s\t%s\n", $seqid, $headt, $tailt;
 	my $curseq = $seq->seq();
 	$seq->seq(('N' x $headt).$curseq.('N' x $tailt));
	$seq_out->write_seq($seq);
}

undef $seq_in;
undef $seq_out;

exit 0;