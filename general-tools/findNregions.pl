#!/usr/bin/perl -w
use strict;
use File::Basename;
use Bio::SeqIO;
use Getopt::Std;

# Search N-regions in multifasta (genomes)
# and produce a BED file with found locations
# adapted from http://stackoverflow.com/questions/10319696
#
# Stephane Plaisance (VIB-NC+BITS) 2015/04/02; v1.0
# visit our Git: https://github.com/BITS-VIB

# disable buffering to get output during long process (loop)
$|=1;

getopts('i:l:h');
our($opt_i, $opt_l, $opt_h);

my $usage="## Usage: findNregions.pl <-i fasta-file>
# Additional optional parameters are:
# <-l minsize in bps (100)>
# <-h to display this help>";

####################
# declare variables
####################

my $fastain = $opt_i || die $usage."\n";
my $minlen = $opt_l || 100;
defined($opt_h) && die $usage."\n";

# open stream from BED file
my $outpath = dirname($fastain);
my $basename = basename($fastain);
(my $outbase = $basename) =~ s/\.[^.]+$//;

# include size limit and max intensity in file names
my $outfile = $outpath."/".$outbase."_N-regions.bed";
open OUT, "> $outfile" || die $!;

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
	print STDERR "## Searching sequence $seqid for $motif\n";
	my $sequence = $seq_obj->seq();

	# scan for motif and report hits
	while ($sequence =~ m/$motif/gi) {
		$counter++;
		my $match_start = $-[0]+1;
		my $match_end = $+[0];
		my $match_seq = $&;
		# print in BED5 format
		print OUT join("\t", $seqid, $match_start,
			$match_end, "N-region", length($&), "+")."\n";
		}

	# end for this sequence
	print STDERR "# found $counter matches for $seqid\n";
	$totcnt += $counter;
	}

print STDERR "#\n# found a total of $totcnt N-regions of $minlen bps or more\n";
close OUT;

exit 0;
