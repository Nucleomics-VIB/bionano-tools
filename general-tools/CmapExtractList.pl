#!/usr/bin/perl

# Extract a list of sequence records (list.txt) from a full.cmap file
# translate sequence names into cmap-IDs
# identify cmap-IDs present in the original cmap
# Save these records to a new .cmap file
#
# Stephane Plaisance (VIB-NC+BITS) 2016/10/20; v1.0
# visit our Git: https://github.com/BITS-VIB

use warnings;
use strict;
use Getopt::Std;
use File::Basename;

my $version = "1.0";

my $usage="## Usage: CmapExtractList.pl 
## script version:".$version."
# <-i cmap_file (required)> 
# <-r record_list file (required)>
# <-k key file (required)>
# <-o outfile | (optional) default from infile name>
# <-h to display this help>";

# disable buffering to get output during long process (loop)
$|=1; 

####################
# declare variables
####################
getopts('i:k:r:o:zh');
our ($opt_i, $opt_k, $opt_r, $opt_o, $opt_h);

# unix | macOS
my $sep = "/";

my $infile = $opt_i || die $usage."\n";
my ($filename,$filedir,$fileext) = fileparse($infile, qr/\.[^.]*/);
my $records = $opt_r || die $usage."\n";
my $keyfile = $opt_k || die $usage."\n";
my ($recname,$recdir,$recext) = fileparse($records, qr/\.[^.]*/);
my $outfile = $opt_o || $filedir.$sep.$recname."-from-".$filename.".cmap";
# simplify if local folder
$outfile =~ s/^\.\/\///;
defined($opt_h) && die $usage."\n";

# load keys from keyfile
print STDOUT "\n# loading record list\n";
open (IDS, $records ) || die "cannot open \"$records\"!";
# load ID-list to array
my @ids = <IDS>;
chomp @ids;
# remove training '_obj'
@ids = map {s/_obj$//g; $_; } @ids;
close IDS;
print STDERR "## ".scalar(@ids)." listed records found\n";

# load keys from keyfile
print STDOUT "\n# loading keypair data\n";
open (KEYS, $keyfile ) || die "cannot open \"$keyfile\"!";
# load keypairs to hash
my %cmapid = ();
while (my $line=<KEYS>) {
	next unless $line =~ /^[0-9]/;
	my @field = split /\t/, $line;
	$cmapid{$field[1]} = $field[0];
}
close KEYS;
my @query = keys %cmapid;
print STDERR "## ".scalar(@query)." listed key pairs found\n";

# find IDs from sequence names
my @filter = ();
my $counter = 0;
foreach (@ids) {
	next unless defined($cmapid{$_});
	$counter++;
	push @filter, $cmapid{$_};
}

# parse cmap file and filter from list
print STDOUT "\n# parsing cmap and filtering data\n";
my %found = ();
open CMAP, $infile or die $!;
open OUT, "> $outfile" || die $!;
while (my $line = <CMAP>) {
	if ($line =~ /^#/) {
		print OUT $line;
		next;
		}
	# this is data
	my @field = split /\t/, $line;
	if ( grep( /^$field[0]$/, @filter ) ) {
	  $found{$field[0]} = 1;
	  print OUT $line;
	}
}
my $cnt = keys %found;
print STDERR "## ".$cnt." cmap records saved to a new file\n";
close CMAP;
close OUT;
exit 0;
