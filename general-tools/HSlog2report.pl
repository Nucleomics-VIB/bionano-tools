#!/usr/bin/perl -w
use strict;

# script: HSlog2report.pl
# convert HYBRID_SCAFFOLD.log to hybrid_scaffold_informatics_report.txt
# the report file is absent when working on IrysSolve and required for BN-access
# @Bionano developpers: this code could be called at the end of the HybridScaffold.pl! 
#
# Stephane Plaisance (VIB-NC+BITS) 2017/05/19; v1.0
# visit our Git: https://github.com/Nucleomics-VIB

# read HYBRID_SCAFFOLD.log in variable, filter, and convert to array
local $/=undef;
open FILE, "HYBRID_SCAFFOLD.log" or die "Couldn't open file: $!";
my $file = <FILE>;
close FILE;

# keep only block of interrest
$file =~ /\n\nCalculating statistics\.\.\.\n\n(Original\ BioNano\ Genome\ Map\ statistics.*)\n\n\nCalculating\ CMAP\ statistics\ complete\ in\ /s;

# convert to array
my @lines = split(/\n/, $1);

# write the converted text to file
open OUT, "> hybrid_scaffold_informatics_report.txt" or die "Couldn't open file: $!";

my $c=0;
foreach (@lines) {
	$c++;
	next if /Running command/;
	next if /^$/;
	print OUT "\n" if ($c gt 1 && /:$/);
	print OUT $_."\n";
}
print OUT "\n";

close OUT;

exit 0;
