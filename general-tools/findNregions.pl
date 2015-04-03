#!/usr/bin/perl -w
use strict;
use File::Basename;
use Bio::SeqIO;
use Getopt::Std;
use File::Tee qw(tee);

# Search N-regions in multifasta (genomes)
# and produce a BED file with found locations
# use the knicker key file to rename contigs
#
# adapted from http://stackoverflow.com/questions/10319696
#
# Stephane Plaisance (VIB-NC+BITS) 2015/04/02; v1.0
# visit our Git: https://github.com/BITS-VIB

# disable buffering to get output during long process (loop)
$|=1;

getopts('i:k:l:h');
our($opt_i, $opt_k, $opt_l, $opt_h);

my $usage="## Usage: findNregions.pl <-i fasta-file> <-k key-file to rename contigs>
# Additional optional parameters are:
# <-l minsize in bps (100)>
# <-h to display this help>";

####################
# declare variables
####################

my $fastain = $opt_i || die $usage."\n";
my $keyfile = $opt_k || die $usage."\n";
my $minlen = $opt_l || 100;
defined($opt_h) && die $usage."\n";

# if key-file was provided, check it and create hash for renaming
our %keyhash = ();
our $present = 0 ;
our $absent = 0 ;

# load renaming table from file
open KEYS, $keyfile or die $!;
while (<KEYS>) {
	chomp;
	next if ! ($_ =~ /^[0-9]/); # ignore header lines
	my ($CompntId, $CompntName, $CompntLength) = split "\t";
	$keyhash{$CompntName} = $CompntId;
}
close KEYS;

# open stream from BED file
my $outpath = dirname($fastain);
my $basename = basename($fastain);
(my $outbase = $basename) =~ s/\.[^.]+$//;

# include size limit and max intensity in file names
my $outfile = $outpath."/".$outbase."-".$minlen."bps_N-regions.bed";
open OUT, "> $outfile" || die $!;

# keep log copy of STDOUT (comment out if you do not have 'File::Tee' installed
tee STDOUT, '>', $outfile."_log.txt" or die $!;

# create parser for multiple fasta files
my $parser = Bio::SeqIO->new(-file => $fastain, -format => 'Fasta');

# look for $minlen N's in a row
my $motif="[N]{".$minlen.",}";
my $totcnt = 0;

############################################
# loop over records and return hits to BED #
############################################

while(my $seq_obj = $parser->next_seq()) {
	my $counter=0;
	# load id and sequence into strings
	my $seqid = $seq_obj->id;

	# check if cmap has this fasta record
	if ( defined $keyhash{$seqid} ) {
		$present += 1;
		print STDOUT "## Searching sequence $seqid for $motif\n";
		my $sequence = $seq_obj->seq();

		# scan for motif and report hits
		while ($sequence =~ m/$motif/gi) {
			$counter++;
			my $match_start = $-[0]+1; # BED is zero-based !
			my $match_end = $+[0];
			my $match_seq = $&;

			# print in BED5 format when present in cmap
			print OUT join("\t", $keyhash{$seqid}, $match_start,
				$match_end, "N-region", length($&), "+")."\n";
			}

		# end for this sequence
		print STDOUT "# found $counter matches for $seqid\n";
		$totcnt += $counter;
		} else {
		 	print STDOUT "# $seqid is absent from the cmap\n";
			$absent += 1;
		}
	}

print STDOUT "#\n# found a total of $totcnt N-regions of $minlen bps or more in $present records\n";
print STDOUT "# $absent records from the original fasta file are absent in the cmap\n";
close OUT;

exit 0;
