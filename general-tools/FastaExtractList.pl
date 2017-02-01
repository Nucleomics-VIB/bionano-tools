#!/usr/bin/perl

# Extract a list of Fasta records (list.txt) from a multifasta file
# Save records to a new Fasta file
#
# Stephane Plaisance (VIB-NC+BITS) 2016/09/22; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Bio::SeqIO;

my $version = "1.0";

my $usage="## Usage: FastaExtractList.pl 
## script version:".$version."
# <-i fasta_file (required)> 
# <-r record_list file (required)>
# <-o outfile | (optional) default from infile name>
# <-z zip results (default OFF)>
# <-h to display this help>";

# disable buffering to get output during long process (loop)
$|=1; 

####################
# declare variables
####################
getopts('i:r:o:zh');
our ($opt_i, $opt_o, $opt_r, $opt_z, $opt_h);

# unix | macOS
my $sep = "/";

my $infile = $opt_i || die $usage."\n";
my ($filename,$filedir,$fileext) = fileparse($infile, qr/\.[^.]*/);
my $records = $opt_r || die $usage."\n";
my ($recname,$recdir,$recext) = fileparse($records, qr/\.[^.]*/);
my $outfile = $opt_o || $filedir.$sep.$recname."-from-".$filename;
# simplify if local folder
$outfile =~ s/^\.\/\///;
my $zipit = defined($opt_z) || undef;
defined($opt_h) && die $usage."\n";

# load keys from keyfile
my @recs = ();

print STDOUT "\n# loading record list\n";
open (IDS, $records ) || die "cannot open \"$records\"!";
# load ID-list to array
my @ids = <IDS>;
chomp @ids;
close IDS;

print STDERR "## ".scalar(@ids)." listed records found\n";

# parse multifasta file and filter from list
my $counter=0;

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

# loop in Fasta and save @ids matching records
while ( my $seq = $seq_in->next_seq() ) {
	#my $first = (split(',', $seq->id))[0];
	my $currid = $seq->id;
	#my @f=split(/\|/, $first);
	#print STDERR $f[1]." : ".$first."\n";

	# search match in @ids
	if (grep(/$currid/, @ids) ) {
		$counter++;
		$counter =~ /00$/ && print STDERR "## ".$counter." records processed\n";
		$seq_out-> write_seq($seq);
		# remove found element to speed up process
		@ids = grep (! /$currid/, @ids);
		#print STDERR scalar(@ids)."\n";
	} 
}

print STDERR "# Found and extracted ".$counter." fasta records matching ".$records." to ".$outfile."\n";

# close fh
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
