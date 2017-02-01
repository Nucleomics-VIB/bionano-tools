#!/usr/bin/perl

# Add N's and nicking site to fasta records shorter that 20k
# desperate attempt to rescue some of the shorter contigs from a NGS assembly
# from discussion with Kees-Jan (THX)
#
# Stephane Plaisance (VIB-NC+BITS) 2016/09/12; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Std;

my $version = "1.0";

my $usage="## Usage: fastaPad20k.pl
## script version:".$version."
# <-i fasta_file (required)>
# <-n nicking motif (required)>
#   Nt.BspQI: GCTCTTCN
#   Nt.BbvCI: CCTCAGC
#   Nb.Bsml: GAATGCN
#   Nb.BbvCI: CCTCAGC
#   Nb.BsrD1: GCAATGNN
#   Nb.BssSI: CACGAG
# <-o outfile | default from infile name>
# <-z zip results (default OFF)>
# <-h to display this help>";

####################
# declare variables
####################
getopts('i:o:n:zh');
our ($opt_i, $opt_o, $opt_n, $opt_z, $opt_h);

my $infile = $opt_i || die $usage."\n";
my $outfile = $opt_o || "20kpadded_".$infile;
my $nickmotif = $opt_n || die $usage."\n";
my $zipit = defined($opt_z) || undef;
defined($opt_h) && die $usage."\n";

# process multifasta
my $seq_in = OpenArchiveFile($infile);
my $seq_out;

if ( defined($zipit) ) {
	my $bgzip = `which bgzip`;
	die "No bgzip command available\n" unless ( $bgzip );
	chomp($bgzip);
	my $fh;
	open $fh,  " | $bgzip -c >  $outfile\.gz" || die $!;
	$seq_out = Bio::SeqIO->new( -format => 'Fasta', -fh => $fh);
} else {
	$seq_out = Bio::SeqIO -> new( -format => 'Fasta', -file => ">$outfile" );
}

while ( my $seq = $seq_in->next_seq() ) {
	# test size or pad
	my $seqlen = $seq->length;
 	if ( $seqlen < 20000 ) {
 		# make it 21k for good measure
 		my $padlen = 21000 - $seqlen;
 		my $padstring = ("N" x $padlen);
 		my $curseq = $seq->seq();
 		$seq->seq($curseq.$padstring.$nickmotif);
 		my $newname = ( $seq->display_id()."_20k-padded" );
 		$seq->display_id($newname);
	}
	$seq_out->write_seq($seq);
}

undef $seq_in;
undef $seq_out;

exit 0;

#### Subs ####
sub OpenArchiveFile {
    my $infile = shift;
    my $FH;
    if ($infile =~ /.fa$|.fasta$|.fna$/i) {
    $FH = Bio::SeqIO -> new(-file => "$infile", -format => 'Fasta');
    }
    elsif ($infile =~ /.fa.bz2$|.fasta.bz2$|.fna.bz2$/i) {
    $FH = Bio::SeqIO -> new(-file => "bgzip -c $infile | ", -format => 'Fasta');
    }
    elsif ($infile =~ /.fa.gz$|.fasta.gz|.fna.gz/i) {
    $FH = Bio::SeqIO -> new(-file => "gzip -cd $infile | ", -format => 'Fasta');
    }
    elsif ($infile =~ /.fa.zip$|.fasta.zip$|.fna.zip$/i) {
    $FH = Bio::SeqIO -> new(-file => "unzip -p $infile | ", -format => 'Fasta');
    } else {
	die ("$!: do not recognise file type $infile");
	# if this happens add, the file type with correct opening proc
    }
    return $FH;
}
