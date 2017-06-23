#!/usr/bin/perl

# Print prettier TEXT & HTML reports from BioNano 'exp_informaticsReport.txt'
# Should take haploid reports as well as diploid ones regardless of the chosen iteration count
# Stephane Plaisance (VIB-NC+BITS) 2017/06/23; v1.1
# visit our Git: https://github.com/BITS-VIB

use warnings;
use strict;
use Getopt::Std;
use File::Basename;

# specials
#use Data::Dumper;
use Array::Transpose;
# @array=transpose(\@array);

my $version ="2017-06-23, version 1.1";

my $usage="## Usage: prettyReport.pl <-i exp_informaticsReport.txt> 
# <-o base-name for output (optional)>
# <-w also create html output (optional)
# <-h to display this help>
# script version: $version";

my $citeus="<i># Stephane Plaisance (VIB-NC) $version\n# visit our Git: <a href=\"https://github.com/Nucleomics-VIB\">Nucleomics-VIB</a></i>";

####################
# declare variables
####################
getopts('i:o:wh');
our ($opt_i, $opt_o, $opt_w, $opt_h);

my $infile = $opt_i || die $usage."\n";
my $outbase = $opt_o || "prettyReport";
my $makehtml = defined($opt_w) || undef;
defined($opt_h) && die $usage."\n";

open FILE, $infile or die $!;
my $outpath = dirname($infile);
my $outfile = $outpath."/".$outbase.".txt";
my $html = $outfile;
$html =~ s/.txt/.html/;

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

if ( defined($makehtml) ) {
	open HTML, "> $html" || die $!;
	html_header("Status Info de-novo Assembly");
}

# print command and version
print OUT "## De Novo Assembly results\n";
print OUT "\n# Assembly Command :\n\'".$command."\'\n";
print OUT "\n# Pipeline version :".$pipelineversion."\n";

if ( defined($makehtml) ) {
	text2title("Status Info de-novo Assembly");
	text2line("HTML version of the <b>exp_informaticsReport.txt</b> file produced during denovo assembly on the Irys Solve server");
	array2lines( split("\n", $citeus) );
	text2subtitle("# Assembly Command");
	array2lines( split("\n", $command) );
	text2line("Pipeline version :$pipelineversion");
	text2line("&nbsp");
}

# print molecule noise parameters & stats results
print OUT "\n# Molecule Stats\n";
print OUT join("\t", "stats", @{@{$stats[0]}[0]}, "cvg (x)") . "\n";
for (my $idx=0; $idx < scalar @stats; $idx++) {
    print OUT join("\n", join("\t", $statst[$idx], @{@{$stats[$idx]}[1]})) . "\n";
}

if ( defined($makehtml) ) {
	text2subtitle("# Molecule Statistics");
	text2line("The first row reports metrics for all molecules submitted to de-novo assembly, 
	the second raw reports sorted molecules and raw coverage considering the provided NGS reference length, 
	while the last row reports molecules and raw coverage considering the optical assembly size.");
	my @table = ();
	push @table, [ "stats", @{@{$stats[0]}[0]}, "cvg (x)" ];
	for (my $idx=0; $idx < scalar @stats; $idx++) {
    	push @table, [ $statst[$idx], @{@{$stats[$idx]}[1]} ];
	}	
	array2table( @table );
	text2line("&nbsp");
}

# print assembly stages results
print OUT "\n# Molecule Statistics\n";
print OUT join("\t", "assembly stages", @{@{$stages[0]}[0]}) . "\n";
for (my $idx=0; $idx < scalar @stages; $idx++) {
    print OUT join("\n", join("\t", $stagest[$idx], @{@{$stages[$idx]}[1]})) . "\n";
}

if ( defined($makehtml) ) {
	text2subtitle("# Assembly Stage Statistics");
	text2line("The following table reports statistics for all assembly stages except the RefineFinal last stages.");
	my @table = ();
	push @table, [ "assembly stages", @{@{$stages[0]}[0]} ];
	for (my $idx=0; $idx < scalar @stages; $idx++) {
    	push @table, [ $stagest[$idx], @{@{$stages[$idx]}[1]} ];
	}	
	array2table( @table );
	text2line("&nbsp");
}

# print final stage results
print OUT "\n# Final Assembly Stage Results\n";
print OUT join("\t", "Final assembly stages", @{@{$final[0]}[0]}) . "\n";
for (my $idx=0; $idx < scalar @final; $idx++) {
    print OUT join("\n", join("\t", $finalt[$idx], @{@{$final[$idx]}[1]})) . "\n";
}

if ( defined($makehtml) ) {
	text2subtitle("# Final Assembly Stage Results");
	text2line("RefineFinal last stages either reported in \'haploid\' mode or in \'haplotype-aware\' mode depending on the selected XML settings.");
	my @table = ();
	push @table, [ "assembly stages", @{@{$final[0]}[0]} ];
	for (my $idx=0; $idx < scalar @final; $idx++) {
    	push @table, [ $finalt[$idx], @{@{$final[$idx]}[1]} ];
	}	
	array2table( @table );
	text2line("&nbsp");
}

# print molecule alignment results
print OUT "\n# Molecule alignments\n";
print OUT join("\t", "alignments", @{@{$aligns[0]}[0]}) . "\n";
for (my $idx=0; $idx < scalar @aligns; $idx++) {
    print OUT join("\n", join("\t", $alignst[$idx], @{@{$aligns[$idx]}[1]})) . "\n";
}

if ( defined($makehtml) ) {
	text2subtitle("# Molecule alignments");
	text2line("Molecule alignment statistics and noise parameters computed against the provide reference or the final optical assembly.");
	my @table = ();
	push @table, [ "alignments", @{@{$aligns[0]}[0]} ];
	for (my $idx=0; $idx < scalar @aligns; $idx++) {
    	push @table, [ $alignst[$idx], @{@{$aligns[$idx]}[1]} ];
	}	
	array2table( @table );
	text2line("&nbsp");
}

# print SV results
print OUT "\n# Structural differences with the provided Reference\n";
# remove outer-spaces in titles
@svres = map {s/^\s+|\s+$//g; $_; } @svres;
for (my $idx=0; $idx < scalar @svres; $idx++) {
    print OUT join("\t", split(":", $svres[$idx])) . "\n";
}
close OUT;

if ( defined($makehtml) ) {
	text2subtitle("# Structural differences with the provided Reference");
	text2line("SV counts as reported by the pipeline.");
	my @table = ();
	for (my $idx=0; $idx < scalar @svres; $idx++) {
		my @row = split(":", $svres[$idx]);
    	push @table, [ @row ];
	}	
	array2table( @table );
	text2line("&nbsp");
}

if ( defined($makehtml) ) {
	html_footer();
}

exit 0;

############### SUBS ################

#################
sub html_header {
my $html_title = shift;
print HTML <<"END_TXT";
<!DOCTYPE html>
<html>
<head>
<title>$html_title</title>
<style>
table, th, td {
    border: 1px solid black;
    border-collapse: collapse;
}
th, td {
    padding: 3px;
}
</style>
</head>
<body>\n
END_TXT
}

#################
sub html_footer {
print HTML <<"END_TXT";
</body>
</html>\n
END_TXT
}

################
sub text2title {
my $txt = shift;
print HTML "<h1>$txt</h1>\n"
}

################
sub text2subtitle {
my $txt = shift;
print HTML "<h2>$txt</h2>\n"
}

################
sub text2line {
my $txt = shift;
print HTML "<p>$txt</p>\n"
}

#################
sub array2lines {
my @a = @_;
print HTML "</p>\n";
foreach my $row (@a) {
	print HTML "$row<br>";
	}
	print HTML "</p>\n";
}

#################
sub array2table {
	my @AoA = @_;
	print HTML "<table>\n";
	foreach my $row (@AoA) {
    	print HTML "\t<tr>\n";
    	foreach my $col (@{$row}) {
    		print HTML "\t\t<td>$col</td>\n";
    	}
    print HTML "\t</tr>\n";
	}
	print HTML "</table>\n";
}
