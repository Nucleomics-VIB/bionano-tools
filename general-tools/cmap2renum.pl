#!/usr/bin/perl -w

# renumber a refefence cmap starting at 1
# create a matching key file
#
# required for cmaps where many fasta records have been filtered out,
# leading to gaps in numbering and cmap id's greater than 100,000
# IrysView hybrid scaffold expects the contigIDs to be numbered up to 100,000
#
# Stephane Plaisance (VIB-NC+BITS) 2015/06/13; v1.00
#
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

# handle command parameters
getopts('c:k:p:h');
our($opt_c, $opt_k, $opt_p, $opt_h);

my $usage = "# Usage: cmap2renum.pl <-c cmap-file> <-k key-file>
# Optional parameters:
# -p <prefix (\'new_\')>
# -h <this help message>\n";

my $cmapfile = $opt_c || die $usage;
my $keyfile = $opt_k || die $usage;
my $prefx = $opt_p || "new_";
defined($opt_h) && die $usage . "\n";

my $outpath = dirname($cmapfile);
my $outbase = basename($cmapfile, ".cmap");
my $outcmap = $prefx.$outbase.".cmap";

my $keybase = basename($keyfile, ".txt");
my $outkey = $prefx.$keybase.".txt";

# result files
open OUT, "> $outkey" || die $!;
open OUT2, "> $outcmap" || die $!;

# load full key data into an array
open KEYS, $keyfile or die $!;
our @orikeys = ();

# store key file data
while (my $line = <KEYS>) {
	# this is header
	if ($line !~ /^[0-9]/) {
		print OUT $line;
		next;
		}

	# store remaining rows in array
	next if $line =~ /^$/;
	push @orikeys, $line;
}

close KEYS;

print STDERR "# loaded ".scalar(@orikeys)." key rows\n";

# load full cmap data into an array
open CMAP, $cmapfile or die $!;
our @oricmaps = ();

# parse data file
while (my $line = <CMAP>) {
	# this is header
	if ($line =~ /^#/) {
		print OUT2 $line;
		next;
		}

	# store remaining rows in array
	last if $line =~ /^$/;
	push @oricmaps, $line;
}

close CMAP;

print STDERR "# loaded ".scalar(@oricmaps)." cmap rows\n";
#################################
# parse both arrays and renumber
our $newid = 1;

foreach my $onekey (@orikeys) {

	# get oriID
	my @keyfield = split("\t", $onekey);
	my $keyid = $keyfield[0];

	# output new key-line
	$keyfield[0] = $newid;
	print OUT join("\t", @keyfield);

	# get matching cmap rows
	my $same = 1;

	while ($same == 1) {
		my $onecmap = shift @oricmaps || last;
		my @cmapfield = split("\t", $onecmap);
		my $cmapid = $cmapfield[0];

		# test
		if ($cmapid eq $keyid) {
			$same = 1;
			$cmapfield[0] = $newid;
			# print new row
			print OUT2 join("\t", @cmapfield);
		} else {
			$same = 0;
			$newid++;
			$cmapfield[0] = $newid;
			# print new row
			print OUT2 join("\t", @cmapfield);
		}
	}
}

# end
