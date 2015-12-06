#!/usr/bin/perl

# Rename multifasta file using BioNano key file
#
# Stephane Plaisance (VIB-NC+BITS) 2015/12/04; v1.0
# supports compressed files (zip, gzip, bgzip)
#
# visit our Git: https://github.com/BITS-VIB

use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Std;

my $usage="## Usage: 
fastaRename.pl <-i fasta_file (required)> <-k key file (required)>
# <-h to display this help>";

####################
# declare variables
####################
getopts('i:o:k:h');
our ($opt_i, $opt_o, $opt_k, $opt_h);

my $infile = $opt_i || die $usage."\n";
my $outfile = $opt_o || "renamed_".$infile;
my $keyfile = $opt_k || die $usage."\n";
defined($opt_h) && die $usage."\n";

# load keys from keyfile
my @header = "";
my @keys = ();
my %translate = ();

print STDOUT "\n# loading key pairs\n";
open KEYS, $keyfile or die $!;
while (my $line = <KEYS>) {
	$line =~ s/\s+$//;
	next if ($line =~ /^#|^$|^CompntId/);
	#next if ($line =~ /^$/);
	#next if ($line =~ /^CompntId/);
	# fill a hash with replacement numbers
	my @keys = split /\t/, $line;
	$translate{$keys[1]} = $keys[0];
	print STDOUT $keys[1]." => ".$translate{$keys[1]}."\n";
}
close KEYS;
print STDOUT "\n";

# process multifasta
my $seq_in = OpenArchiveFile($infile);
my $seq_out = Bio::SeqIO -> new( -format => 'Fasta', -file => ">$outfile" );
# my $seq_out = Bio::SeqIO -> new( -format => 'Fasta', -file => " | gzip -c >$outfile" );

while ( my $seq = $seq_in->next_seq() ) {
	my $curname = $seq->display_id()." ".$seq->desc;
	my $newname = $translate{$curname};
	print STDOUT "# renaming: \"".$curname."\" to ".$newname."\n";
	$seq->display_id($newname);
	$seq->accession_number("");
	$seq->desc("");
	$seq_out->write_seq($seq); 
	}

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
