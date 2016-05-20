#!/usr/bin/perl -w

# smap2bed.pl (first version: 2016)
# convert the ref/query part of a smap (v0.4) to BED5 (strand is set to '.')
# default use confidence value (9) for the score
# use aligned part of the reference for coordinates
#
# Stephane Plaisance (VIB-NC+BITS) 2016/05/20; v1.0
# initial version
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Statistics::Descriptive;
use List::Util qw( min max );
use List::MoreUtils qw( uniq );

# handle command parameters
getopts('i:x:c:n:s:p:h');
our($opt_i, $opt_x, $opt_c, $opt_n, $opt_s, $opt_p, $opt_h);

my $usage = "Aim: Convert smap data to BED5. You must provide a smap file with -i
# Usage: smap2bed.pl <-i smap-file>
# Optional parameters (smap v0.4) :
# -x <minimal value for Confidence score (default=-1)>
# -c <coordinate system used <'q'=query/'r'=ref> (default='r')
# -n <field number for BED-name (1-based; default to SmapEntryID=1)>
#        1:SmapEntryID 2:QryContigID 3:RefcontigID1 4:RefcontigID2 5:QryStartPos 6:QryEndPos
#        7:RefStartPos 8:RefEndPos 9:Confidence 10:Type 11:XmapID1 12:XmapID2 13:LinkID
#       14:QryStartIdx 15:QryEndIdx 16:RefStartIdx 17:RefEndIdx
# -s <field number for BED-score (1-based; default to Confidence=9)>
#        1:SmapEntryID 2:QryContigID 3:RefcontigID1 4:RefcontigID2 5:QryStartPos 6:QryEndPos
#        7:RefStartPos 8:RefEndPos 9:Confidence 10:Type 11:XmapID1 12:XmapID2 13:LinkID
#       14:QryStartIdx 15:QryEndIdx 16:RefStartIdx 17:RefEndIdx
# -p <percentile for Confidence distribution (default=95>)
# <-h to display this help>";

# // tests for defined-ness rather than truth
# allows user inputing 0 as value where || would not do
# http://stackoverflow.com/questions/1609060/how-can-i-set-default-values-using-getoptstd
defined($opt_h) && die $usage . "\n";
my $inputfile = $opt_i || die $usage;
my $minscore = $opt_x // -1; # '-1' is used for a number of SV types
my $coordinate = $opt_c || "r";
my $namefield = $opt_n || 1;
my $scorefield = $opt_s || 9;
my $percentile = $opt_p // 95;

# test input
grep( /^$coordinate$/, ( "q","r" ) ) || die "-c should be of 'q'/'r'\n";
grep( /^$namefield$/, ( 1..17 ) ) || die "-n should be in [1..17]\n";
grep( /^$scorefield$/, ( 1..17 ) ) || die "-s should be in [1..17]\n";
($minscore >= -1) || die "-x should be in [-1..]\n";
($percentile >= 0 && $percentile <= 100) || die "-p should be in [0..100]\n";

our %fieldnames = (
	1 => "SmapEntryID",
	2 => "QryContigID",
	3 => "RefcontigID1",
	4 => "RefcontigID2",
	5 => "QryStartPos",
	6 => "QryEndPos",
	7 => "RefStartPos",
	8 => "RefEndPos",
	9 => "Confidence",
	10 => "Type",
	11 => "XmapID1",
	12 => "XmapID2",
	13 => "LinkID",
	14 => "QryStartIdx",
	15 => "QryEndIdx",
	16 => "RefStartIdx",
	17 => "RefEndIdx"
	);

# report choices
print STDOUT "\n##### BED-field options #####\n";
print STDOUT "| coordinates: ".$coordinate."\n";
print STDOUT "| seqlab: ".($coordinate eq 'q' ? $fieldnames{2} : ($fieldnames{3}.",".$fieldnames{4}))."\n";
print STDOUT "| start: ".($coordinate eq 'q' ? $fieldnames{5} : $fieldnames{7})."\n";
print STDOUT "| end: ".($coordinate eq 'q' ? $fieldnames{6} : $fieldnames{8})."\n";
print STDOUT "| name: ".$fieldnames{$namefield}."\n";
print STDOUT "| score: ".$fieldnames{$scorefield}."\n";
print STDOUT "| strand: "."."."\n";
print STDOUT "\n# filtering parameters\n";
print STDOUT "| cutoff score: ".$minscore."\n";
print STDOUT "| percentile score: ".$percentile."\n";

# load smap header and process content
open FILE, $inputfile or die $!;
my $outpath = dirname($inputfile);
my $outbase = basename($inputfile, ".smap");
my $outfile = $outbase."_gt".$minscore."_".$coordinate.".bed";

# result files
open OUT, "> $outfile" || die $!;

# declare variables
my $countsmap = 0;
our @comments = ();
our @header = ();
our @colnames = ();
our @coltypes = ();
our @Confidences_all = ();
our @Confidences_filt = ();
# count per type in a hash
our %type_all = ();
our %type_filt = ();

# parse data file
while (my $line = <FILE>) {
	if ($line =~ /^#/) {
		parseheader($line);
		next;
		}

	# this is data
	$countsmap++;

	my @field = ( undef, (split /\t/, $line) );
	# further decompose
	my @seqlab = ($coordinate eq 'q' ? $field[2] : ( uniq ( $field[3], $field[4] ) ));
	# ignore $seqlab == '-1' for partial calls
	@seqlab = grep { $_ != '-1' } @seqlab;
	my $start = ($coordinate eq 'q' ? $field[5] : $field[7]);
	my $end = ($coordinate eq 'q' ? $field[6] : $field[8]);
	my $coordstart = int( min($start, $end)+0.5);
	my $coordend = int( max($start, $end)+0.5);
	my $score = $field[$scorefield];

	# print next rows until contig ends
	# all records
	push @Confidences_all, $score;
	$type_all{$field[10]}++;

	# test $field[$scorefield]>$minscore and ?print
	if ($score >= $minscore) {
		# store this score for stats
		push @Confidences_filt, $score;
		$type_filt{$field[10]}++;

		# output data
		my $onerec;
		foreach $onerec (@seqlab) { 
			print OUT join("\t", 
				$onerec,
				$coordstart,
				$coordend,
				$field[$namefield],
				$field[$scorefield],
				"."
				)."\n";
		}
	}
}

close FILE;
close OUT;

######################
### report findings

# get stats for the Confidence scores above threshold
my @quantheader=("min","25%","median","75%","max","mean", $percentile."%");
#report_stats("Confidence scores for All SV calls", \@Confidences_all);
#report_stats("Confidence scores for SV calls above ".$minscore, \@Confidences_filt);

my @all_scores = get_stats(@Confidences_all);
my @filt_scores = get_stats(@Confidences_filt);

# print summary
print STDOUT "\n##### smap header information #####\n";
print STDOUT "| input: ".basename($inputfile)."\n";
print STDOUT "| records: ".$countsmap."\n";
print STDOUT "| Colnames: ";
print STDOUT join(", ", @colnames)."\n";

print STDOUT "\n##### Headers #####\n";
print STDOUT "|".join("\n|", @header)."\n";

print STDOUT "\n##### Comments #####\n";
print STDOUT "|".join("\n|", @comments)."\n";
print STDOUT "##############################\n";

print STDOUT "\n##### Confidence score distribution #####\n";
print STDOUT "# For All calls\n";
print STDOUT "|".join("\t", @quantheader)."\n";
print STDOUT "|".join("\t", @all_scores)."\n";
print STDOUT "# For Calls above cutoff (".$minscore.")\n";
print STDOUT "|".join("\t", @quantheader)."\n";
print STDOUT "|".join("\t", @filt_scores)."\n";
print STDOUT "##############################\n";

print STDOUT "\n##### SV call types #####\n";
print STDOUT "### For All calls\n";
foreach my $str (sort keys %type_all) {
	printf STDOUT "|%-31s %s\n", $str, $type_all{$str};
	}
print STDOUT "### For Calls with Confidence above cutoff (".$minscore.")\n";
foreach my $str (sort keys %type_filt) {
	printf STDOUT "|%-31s %s\n", $str, $type_filt{$str};
	}

################
##### Subs #####

################
sub parseheader {
my $line = shift;
chomp($line);
$line =~ s/\.\.\.$//;
# put header in array
my @arr = split("\t", $line);
if ($arr[0] =~ /^#h/) {
		@colnames = @arr[1 .. $#arr];
	} elsif ($arr[0] =~ /^#f/) {
		@coltypes = @arr[1 .. $#arr];
	} elsif ($#arr == 0) {
		push @comments, $line;
	} else {
		push @header, $line;
	}
}

##############
sub get_stats {
# uses https://metacpan.org/module/Statistics::Descriptive
# test input
return undef unless ( scalar(@_) );

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@_);
my @result = ( 
	sprintf("%.2f", $stat->quantile(0)),
	sprintf("%.2f",$stat->quantile(1)),
	sprintf("%.2f",$stat->quantile(2)),
	sprintf("%.2f",$stat->quantile(3)),
	sprintf("%.2f",$stat->quantile(4)),
	sprintf("%.2f",$stat->mean()),
	sprintf("%.2f",scalar($stat->percentile($percentile)))
	);

return @result;
}

