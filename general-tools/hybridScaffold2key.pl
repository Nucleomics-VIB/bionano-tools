#!/usr/bin/perl

# Take <NGScontigs_HYBRID_SCAFFOLD.fasta> file and associated <NGScontigs_HYBRID_SCAFFOLD.xmap_sorted.xmap>
# compute sequence lengths from Fasta
# extract BioNano ID from xmap
#  (in same appearance order as Fasta)
# create a BioNano key file for name to ID translation 
# REM: normally, this script is useless since the key is the numeric part XXX at the end of each SuperScaffold_XXX name
#
# Stephane Plaisance (VIB-NC+BITS) 2017/06/27; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Bio::SeqIO;

my $version = "1.0, 2017-06-27";

my $usage="## Usage: hybridScaffold2key.pl 
# <-i NGScontigs_HYBRID_SCAFFOLD.fasta (required)> 
# <-x NGScontigs_HYBRID_SCAFFOLD.xmap_sorted.xmap (required)>
# <-h to display this help>
## script version:".$version;

####################
# declare variables
####################
getopts('i:x:h');
our ($opt_i, $opt_x, $opt_h);

# unix | macOS
my $sep = "/";

my $fasta = $opt_i || die $usage."\n";
my ($fastaname,$fastadir,$fastaext) = fileparse($fasta, qr/\.[^.]*/);
my $xmap = $opt_x || die $usage."\n";
my ($xname,$xdir,$xext) = fileparse($xmap, qr/\.[^.]*/);
my $keyfile = $fastadir.$sep.$fastaname."_keys.txt";
# simplify if local folder
$keyfile =~ s/^\.\/\///;
defined($opt_h) && die $usage."\n";

# open output handle
open KEY, "> $keyfile" || die $!;
printf KEY "%s\t%s\t%s\n", "CompntId", "CompntName", "CompntLength";

print STDOUT "\n# loading ID list from xmap\n";
open XMAP, $xmap or die $!;
my @all = ();

while (my $line = <XMAP>) {
	# ignore comment header lines
	next if $line =~ /^#/;
	my @field = split "\t", $line;
	my $id = $field[2];
	push @all, $id;
}
close XMAP;

# keep unique IDs
my @ids = uniques(@all);
print STDERR "## ".scalar(@ids)." found\n";

# parse multifasta file and filter from list
my $seq_in = OpenArchiveFile($fasta);
my $seqnum = 0;

# loop in Fasta and save @ids matching records
while ( my $seq = $seq_in->next_seq() ) {
	my $currname = $seq->id;
	my $currlen = $seq->length;
	my $xmapID = $ids[$seqnum];
	printf KEY "%s\t%s\t%s\n", $xmapID, $currname, $currlen;
	$seqnum++;
}

# close fh
undef $seq_in;
print STDERR "# Found and processed ".$seqnum." fasta records to ".$keyfile."\n";

exit 0;

#### Subs ####
sub uniques {
my @data = @_;
my @unique;
my %seen;
foreach my $value (@data) {
  if (! $seen{$value}++ ) {
    push @unique, $value;
  	}
}
return @unique;
}

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
