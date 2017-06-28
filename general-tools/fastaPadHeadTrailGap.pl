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

my $version = "1.0";

my $usage="## Usage: fastaPadHeadTrailGap.pl
## script version:".$version."
# <-i NGScontigs_HYBRID_SCAFFOLD.fasta file (required)>
# <-t NGScontigs_HYBRID_SCAFFOLD_trimHeadTailGap.coord file (required)>
# <-h to display this help>";

####################
# declare variables
####################
getopts('i:t:h');
our ($opt_i, $opt_t, $opt_h);

my $infile = $opt_i || die $usage."\n";
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

###########################################
# process NGScontigs_HYBRID_SCAFFOLD.fasta 
###########################################

my $seq_in = Bio::SeqIO -> new(-file => "$infile", -format => 'Fasta');
my $seq_out = Bio::SeqIO -> new(-file => "> untrimmed_$infile", -format => 'Fasta');

while ( my $seq = $seq_in->next_seq() ) {
	my $seqid = $seq->display_id();
	print STDOUT "# processing \'$seqid\' : ";
	print STDOUT "adding \'$headtrim{$seqid}\' heading Ns, ";
	print STDOUT "and \'$tailtrim{$seqid}\' trailing Ns\n";
 	my $headt = 'N' x $headtrim{$seqid};
 	my $tailt = 'N' x $tailtrim{$seqid};
 	my $curseq = $seq->seq();
 	$seq->seq($headt.$curseq.$tailt);
	$seq_out->write_seq($seq);
}

undef $seq_in;
undef $seq_out;

exit 0;