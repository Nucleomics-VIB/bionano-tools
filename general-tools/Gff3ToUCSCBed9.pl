#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

# convert GFF3 annotations to UCSC bed-9 format for Bionano access
# also convert original contigs names to BioNano IDs from a matching key file
#
# Stephane Plaisance (VIB-NC) initial version 2017-06-16; v1.0
# visit our Git: https://github.com/Nucleomics-VIB

my $version ="1.0, 2017-06-16";

my $usage="# Gff3ToUCSCBed9.pl (version $version)
## Usage: Gff3ToUCSCBed9.pl <-i gff3 file (required)> <-k key file (required)>
# optional <-o name for bed output>
# optional <-f feature type (default to gene, can be exon, cds, ...)>
# <-h to display this help>";

####################
# declare variables
####################
getopts('i:o:k:f:h');
our ($opt_i, $opt_o, $opt_k, $opt_f, $opt_h);

my $infile = $opt_i || die $usage."\n";
my $keyfile = $opt_k || die $usage."\n";
my $feature = $opt_f || "gene";
my @extensions = ( ".gff", ".gff3", ".gff.gz", ".gff3.gz", ".gff.zip", ".gff3.zip", ".gff.bz2", ".gff3.bz2" );
my $outfile = $opt_o || basename($infile, @extensions)."-".$feature.".bed";
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
print STDOUT "\n";

close KEYS;

# process GFF annotation file
my $gff_in = OpenArchiveFile($infile);
open OUT, "> $outfile" || die $!;
my $cnt = 0;

while ( my $line = <$gff_in> ) {
	# pass comment lines
	next if $line =~ /^#/;
	# split in array
	my @data = split("\t", $line);
	# is of feature type & contig/chr has a key?
	if ( $data[2] =~ m/$feature/i && defined($translate{$data[0]}) ) {
		# replace by key ID
		$data[0] = $translate{$data[0]};
		$cnt++;
		#print STDOUT join("\t", @data)."\n";
		convert_to_bed($cnt, @data);
		}
	}
print STDOUT "# ".$cnt." features of type \'".$feature."\' have been converted\n";
print STDOUT "# consider sorting your results with \'cat $outfile | sort -k1n,1 -k2n,2 > sorted_$outfile\'\n";

undef $gff_in;

exit 0;

#### Subs ####
sub convert_to_bed {
my ($cnt, @data) = @_;
my @bed = ();
my $name = "na";
# extract name from col9
if ($data[8] =~ /(^.*;Name=)([^;]*)(;.*$)/) {
	$name = $2;
	}
# order BED fields
$bed[0] = $data[0];
$bed[1] = $data[3]-1;
$bed[2] = $data[4];
$bed[3] = $name;
$bed[4] = $cnt;
$bed[5] = $data[6];
$bed[6] = $data[3]-1;
$bed[7] = $data[4];
$bed[8] = "10,0,255";
# output results
print OUT join("\t", @bed)."\n";
}

sub OpenArchiveFile {
    my $infile = shift;
    my $FH;
    if ($infile =~ /\.gff(.)?$/i) {
    	open ($FH, $infile) or die $!;
	    }
    elsif ($infile =~ /\.gff(.)?.bz2$/i) {
    	open ($FH, "bgzip -c $infile | ") or die $!;
    	}
    elsif ($infile =~ /\.gff(.)?.gz/i) {
    	open ($FH, "gzip -cd $infile | ") or die $!;
    	}
    elsif ($infile =~ /\.gff(.)?.zip$/i) {
    	open ($FH, "unzip -p $infile | ")  or die $!;
    } else {
		die ("$!: do not recognise GFF file extension in $infile");
		# if this happens add, the file type with correct opening proc
    }
    return $FH;
}

#### example
# ensembl GFF3 for hg19
# ##gff-version   3
# ##sequence-region   1 1 249250621
# ## <more such lines>
# 
# 1	GRCh37	chromosome	1	249250621	.	.	.	ID=chromosome:1;Alias=CM000663.1,NC_000001.10
# ###
# 1	cpg	biological_region	10469	11240	1.3e+03	.	.	external_name=oe %3D 0.79;logic_name=cpg
# 1	Eponine	biological_region	10650	10657	0.999	+	.	logic_name=eponine
# 1	Eponine	biological_region	10656	10658	0.999	-	.	logic_name=eponine
# 1	Eponine	biological_region	10678	10687	0.999	+	.	logic_name=eponine
# 1	Eponine	biological_region	10682	10689	0.999	-	.	logic_name=eponine
# 1	Eponine	biological_region	10707	10716	0.999	+	.	logic_name=eponine
# 1	Eponine	biological_region	10709	10719	0.999	-	.	logic_name=eponine
# 1	Eponine	biological_region	10736	10748	0.999	-	.	logic_name=eponine
# 1	Eponine	biological_region	10737	10744	0.999	+	.	logic_name=eponine
# 1	Eponine	biological_region	10766	10773	0.999	+	.	logic_name=eponine
# 1	Eponine	biological_region	10771	10780	0.999	-	.	logic_name=eponine
# 1	Eponine	biological_region	10796	10801	0.999	+	.	logic_name=eponine
# 1	Eponine	biological_region	10811	10820	0.999	-	.	logic_name=eponine
# 1	Eponine	biological_region	10870	10872	0.999	+	.	logic_name=eponine
# 1	Eponine	biological_region	10890	10894	0.999	-	.	logic_name=eponine
# 1	ensembl_havana	pseudogene	11869	14412	.	+	.	ID=gene:ENSG00000223972;Name=DDX11L1;biotype=pseudogene;description=DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1 [Source:HGNC Symbol%3BAcc:37102];gene_id=ENSG00000223972;logic_name=ensembl_havana_gene;version=4
# 1	havana	processed_transcript	11869	14409	.	+	.	ID=transcript:ENST00000456328;Parent=gene:ENSG00000223972;Name=DDX11L1-002;biotype=processed_transcript;havana_transcript=OTTHUMT00000362751;havana_version=1;tag=basic;transcript_id=ENST00000456328;version=2
# 1	havana	exon	11869	12227	.	+	.	Parent=transcript:ENST00000456328;Name=ENSE00002234944;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00002234944;rank=1;version=1
# 1	havana	exon	12613	12721	.	+	.	Parent=transcript:ENST00000456328;Name=ENSE00003582793;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003582793;rank=2;version=1
# 1	havana	exon	13221	14409	.	+	.	Parent=transcript:ENST00000456328;Name=ENSE00002312635;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00002312635;rank=3;version=1
# 1	ensembl	pseudogenic_transcript	11872	14412	.	+	.	ID=transcript:ENST00000515242;Parent=gene:ENSG00000223972;Name=DDX11L1-201;biotype=transcribed_unprocessed_pseudogene;transcript_id=ENST00000515242;version=2
# 1	ensembl	exon	11872	12227	.	+	.	Parent=transcript:ENST00000515242;Name=ENSE00002234632;constitutive=0;ensembl_end_phase=2;ensembl_phase=-1;exon_id=ENSE00002234632;rank=1;version=1
# 1	ensembl	exon	12613	12721	.	+	.	Parent=transcript:ENST00000515242;Name=ENSE00003608237;constitutive=0;ensembl_end_phase=0;ensembl_phase=2;exon_id=ENSE00003608237;rank=2;version=1
# 1	ensembl	exon	13225	14412	.	+	.	Parent=transcript:ENST00000515242;Name=ENSE00002306041;constitutive=0;ensembl_end_phase=-1;ensembl_phase=0;exon_id=ENSE00002306041;rank=3;version=1
# 1	ensembl	pseudogenic_transcript	11874	14409	.	+	.	ID=transcript:ENST00000518655;Parent=gene:ENSG00000223972;Name=DDX11L1-202;biotype=transcribed_unprocessed_pseudogene;transcript_id=ENST00000518655;version=2
# 1	ensembl	exon	11874	12227	.	+	.	Parent=transcript:ENST00000518655;Name=ENSE00002269724;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00002269724;rank=1;version=1
# 1	ensembl	exon	12595	12721	.	+	.	Parent=transcript:ENST00000518655;Name=ENSE00002270865;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00002270865;rank=2;version=1
# 1	ensembl	exon	13403	13655	.	+	.	Parent=transcript:ENST00000518655;Name=ENSE00002216795;constitutive=0;ensembl_end_phase=0;ensembl_phase=-1;exon_id=ENSE00002216795;rank=3;version=1
# 1	ensembl	exon	13661	14409	.	+	.	Parent=transcript:ENST00000518655;Name=ENSE00002303382;constitutive=0;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSE00002303382;rank=4;version=1
# 1	havana	pseudogenic_transcript	12010	13670	.	+	.	ID=transcript:ENST00000450305;Parent=gene:ENSG00000223972;Name=DDX11L1-001;biotype=transcribed_unprocessed_pseudogene;havana_transcript=OTTHUMT00000002844;havana_version=2;transcript_id=ENST00000450305;version=2
# 1	havana	exon	12010	12057	.	+	.	Parent=transcript:ENST00000450305;Name=ENSE00001948541;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001948541;rank=1;version=1
# 1	havana	exon	12179	12227	.	+	.	Parent=transcript:ENST00000450305;Name=ENSE00001671638;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001671638;rank=2;version=2
# 1	havana	exon	12613	12697	.	+	.	Parent=transcript:ENST00000450305;Name=ENSE00001758273;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001758273;rank=3;version=2
# 1	havana	exon	12975	13052	.	+	.	Parent=transcript:ENST00000450305;Name=ENSE00001799933;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001799933;rank=4;version=2
# 1	havana	exon	13221	13374	.	+	.	Parent=transcript:ENST00000450305;Name=ENSE00001746346;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001746346;rank=5;version=2
# 1	havana	exon	13453	13670	.	+	.	Parent=transcript:ENST00000450305;Name=ENSE00001863096;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001863096;rank=6;version=1
# ###

# # Bionano hg19 annotations
# 1	11873	14409	DDX11L1	1	+	11873	14409	10,0,255
# 1	14361	19759	WASH7P	2	+	14361	19759	10,0,255
# 1	14406	29370	WASH7P	3	+	14406	29370	10,0,255
# 1	34610	36081	FAM138F	4	+	34610	36081	10,0,255
# 1	69090	70008	OR4F5	5	+	69090	70008	10,0,255
# 1	134772	140566	LOC729737	6	+	134772	140566	10,0,255
# 1	321083	321115	DQ597235	7	+	321083	321115	10,0,255
# 1	321145	321207	DQ599768	8	+	321145	321207	10,0,255
# 1	322036	326938	LOC100133331	9	+	322036	326938	10,0,255
# 1	327545	328439	LOC388312	10	+	327545	328439	10,0,255
# 1	367658	368597	OR4F29	11	+	367658	368597	10,0,255
# 1	420205	421839	BC036251	12	+	420205	421839	10,0,255
# 1	566092	566115	JA429830	13	+	566092	566115	10,0,255

