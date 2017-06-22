#!/usr/bin/perl

# Print prettier TEXT report from BioNano 'exp_informaticsReport.txt'
# Should take haploid reports as well as diploid ones regardless fo the chosen iteration count
# Stephane Plaisance (VIB-NC+BITS) 2017/06/20; v1.0
# visit our Git: https://github.com/BITS-VIB

use warnings;
use strict;
use Getopt::Std;
use File::Basename;

# specials
#use Data::Dumper;
use Array::Transpose;
# @array=transpose(\@array);

my $usage="## Usage: 
prettyReport.pl <-i exp_informaticsReport.txt> <-o base-name for output (optional)>
# <-h to display this help>";

####################
# declare variables
####################
getopts('i:o:h');
our ($opt_i, $opt_o, $opt_h);

my $infile = $opt_i || die $usage."\n";

my $outbase = $opt_o || "prettyReport";
defined($opt_h) && die $usage."\n";

open FILE, $infile or die $!;
my $outpath = dirname($infile);
my $outfile = $outpath."/".$outbase.".txt";

# declare variables
my $command;
my $pipelineversion;
my @stats;
my @statst;
my @stages;
my @stagest;
my @final;
my @finalt;
my @aligns;
my @alignst;
my @svres;

while (my $line = <FILE>) {

	if ($line =~ m/pipelineCL.py/) {
		chomp($line);
		$command = $line;
		$command =~ s/\ -/\n\ -/g;
		readline(FILE);
		}
	
	if ($line =~ m/Informatics\ Report\ Version\:/) {
		chomp($line);
		$pipelineversion = (split "\:", $line)[1];
		readline(FILE);
		}

	if ($line =~ m/Reading\ molecule\ stats\ from\ (.*)\:$/) {
		chomp($line);
		my $source = (split "/", $1)[-1];
		push @statst, $source;
		readline(FILE);
		# Reading molecule stats from /home/bionano/data/2312_ArabicaWS-2.5.1/2312-Arabica_FM/output/all.bnx:
		# read until empty line then store
		my @data = ();

		while ( $line = <FILE> ) {
			chomp($line);
			my ($title, $value);
			last if $line =~ /^$/;
			my @fields =  split(":", $line);
			($title = $fields[0]) =~ s/^\s+|\s+$//g;
			($value = $fields[1]) =~ s/^\s+|\s+$//g;
			push @data, [( $title, $value )];
			}
		@data = transpose(\@data);
		push @stats, [@data];
		}

	if ($line =~ m/Stage\ Summary\:\ Characterize(.*\ refineFinal1)$/) {
		# refineFinal1 stages
		chomp($line);
		my $source = $1;
		push @finalt, $source;
		# Stage Summary: CharacterizeDefault|Final refineFinal1
		# read until empty line then store
		my @data = ();

		while ( $line = <FILE> ) {
			chomp($line);
			my ($title, $value);
			last if $line =~ /^$/;
			my @fields =  split(":", $line);
			($title = $fields[0]) =~ s/^\s+|\s+$//g;
			($value = $fields[1]) =~ s/^\s+|\s+$//g;
			push @data, [( $title, $value )];
			}
		@data = transpose(\@data);
		push @final, [@data];
	} 
	
	if ($line =~ m/Stage\ Summary\:\ CharacterizeDefault\ (.*)$/) {
		# other stages
		chomp($line);
		my $source = $1;
		push @stagest, $source;
		# Stage Summary: CharacterizeDefault XXXXX
		# read until empty line then store
		my @data = ();

		while ( $line = <FILE> ) {
			chomp($line);
			my ($title, $value);
			last if $line =~ /^$/;
			my @fields =  split(":", $line);
			($title = $fields[0]) =~ s/^\s+|\s+$//g;
			($value = $fields[1]) =~ s/^\s+|\s+$//g;
			push @data, [( $title, $value )];
			}
		@data = transpose(\@data);
		push @stages, [@data];
		}

	if ($line =~ m/Molecules\ Aligned\ to\ (.*)\:/) {
		# molecule alignment tests
		chomp($line);
		my $source = $1;
		push @alignst, $source;
		# Alignment summaries
		# read until empty line then store
		my @data = ();

		while ( $line = <FILE> ) {
			chomp($line);
			my ($title, $value);
			last if $line =~ /^$/;
			my @fields =  split(":", $line);
			($title = $fields[0]) =~ s/^\s+|\s+$//g;
			($value = $fields[1]) =~ s/^\s+|\s+$//g;
			push @data, [( $title, $value )];
			}
		@data = transpose(\@data);
		push @aligns, [@data];
		}

	if ($line =~ m/SV\ detect\:\ svdetect_exp_refineFinal1_sv/) {
		# SV block as-is
		readline(FILE);
		my @data = ();
		while ( $line = <FILE> ) {
		chomp($line);
		last if $line =~ /^$/;
		push @data, $line;
		}
		push @svres, @data;
	}
}

#print Dumper \@statst;
#print Dumper \@stats;
#print Dumper \@stagest;
#print Dumper \@stages;
#print Dumper \@alignst;
#print Dumper \@aligns;
#print Dumper \@svres;

close FILE;

########################
# print results to OUT #
########################

open OUT, "> $outfile" || die $!;

# print command and version
print OUT "## De Novo Assembly results\n";
print OUT "\n# Assembly Command :\n\'".$command."\'\n";
print OUT "\n# Pipeline version :".$pipelineversion."\n";


# print molecule noise parameters & stats results
print OUT "\n# Molecule Stats\n";
print OUT join("\t", "stats", @{@{$stats[0]}[0]}, "cvg (x)") . "\n";
for (my $idx=0; $idx < scalar @stats; $idx++) {
    print OUT join("\n", join("\t", $statst[$idx], @{@{$stats[$idx]}[1]})) . "\n";
}

# print assembly stages results
print OUT "\n# Assembly Stage Results\n";
print OUT join("\t", "assembly stages", @{@{$stages[0]}[0]}) . "\n";
for (my $idx=0; $idx < scalar @stages; $idx++) {
    print OUT join("\n", join("\t", $stagest[$idx], @{@{$stages[$idx]}[1]})) . "\n";
}

# print final stage results
print OUT "\n# Final Assembly Stage Results\n";
print OUT join("\t", "Final assembly stages", @{@{$final[0]}[0]}) . "\n";
for (my $idx=0; $idx < scalar @final; $idx++) {
    print OUT join("\n", join("\t", $finalt[$idx], @{@{$final[$idx]}[1]})) . "\n";
}

# print molecule alignment results
print OUT "\n# Molecule alignments\n";
print OUT join("\t", "alignments", @{@{$aligns[0]}[0]}) . "\n";
for (my $idx=0; $idx < scalar @aligns; $idx++) {
    print OUT join("\n", join("\t", $alignst[$idx], @{@{$aligns[$idx]}[1]})) . "\n";
}

# print SV results
print OUT "\n# Structural differences with the provided Reference\n";
# remove outer-spaces in titles
@svres = map {s/^\s+|\s+$//g; $_; } @svres;
for (my $idx=0; $idx < scalar @svres; $idx++) {
    print OUT join("\t", split(":", $svres[$idx])) . "\n";
}
close OUT;

exit 0;
