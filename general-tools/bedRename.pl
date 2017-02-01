#!/usr/bin/perl

# Rename BED chromosomes using BioNano key file
# only rows with chromosomes present in the key file are converted and saved
#
# Stephane Plaisance (VIB-NC+BITS) 2015/12/04; v1.0
# supports compressed files (zip, gzip, bgzip)
#
# visit our Git: https://github.com/Nucleomics-VIB

use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Std;

my $usage="## Usage: 
bedRename.pl <-i bed file (required)> <-k key file (required)>
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

# disable buffering to get output during long process (loop)
$|=1; 

# load keys from keyfile
my @header = "";
my @keys = ();
my %translate = ();

print STDOUT "\n# loading key pairs\n";
open KEYS, $keyfile or die $!;
while (my $line = <KEYS>) {
	$line =~ s/\s+$//;
	next if ($line =~ /^#|^$|^CompntId/);
	# fill a hash with replacement numbers
	my @keys = split /\t/, $line;
	$translate{$keys[1]} = $keys[0];
	print STDOUT $keys[1]." => ".$translate{$keys[1]}."\n";
}
close KEYS;
print STDOUT "\n";

# process multifasta
my $bed_in = OpenArchiveFile($infile);
open OUT, "> $outfile" || die $!;

while ( my $line = <$bed_in> ) {
	my @fields = split("\t", $line);
	my $cur = $fields[0];
	#print STDOUT $cur," ",(defined($translate{$cur})?$translate{$cur}:"unef")."\n";
	if ( defined($translate{$cur}) ) {
		$fields[0] = $translate{$cur};
		print OUT join("\t", @fields);
		}
	}

close $bed_in;

exit 0;

#### Subs ####
sub OpenArchiveFile {
    my $infile = shift;
    my $FH;
    if ($infile =~ /.bed$/i) {
    	open ($FH, $infile) or die $!;
	    }
    elsif ($infile =~ /.bed.bz2$/i) {
    	open ($FH, "bgzip -c $infile | ") or die $!;
    	}
    elsif ($infile =~ /.bed.gz/i) {
    	open ($FH, "gzip -cd $infile | ") or die $!;
    	}
    elsif ($infile =~ /.bed.zip$/i) {
    	open ($FH, "unzip -p $infile | ")  or die $!;
    } else {
		die ("$!: do not recognise BED file extension in $infile");
		# if this happens add, the file type with correct opening proc
    }
    return $FH;
}
