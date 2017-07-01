#!/usr/bin/perl -w

## dovetail2agp.pl
## first version: 2017-07-01
## extract original coordinate from Dovetail .table.txt
## add header annotations (optional)
## add GAP rows between placed contigs
## save in AGP format

# Stephane Plaisance (VIB-NC+BITS) 2017/01/01; v1.0
# visit our Git: https://github.com/BITS-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Sort::Key::Natural qw(natsort);
#use Data::Dumper;

############################
# handle command parameters
############################

getopts('i:o:n:d:s:D:h');
our ( $opt_i, $opt_o, $opt_n, $opt_d, $opt_s, $opt_D, $opt_h );

my $usage = "Aim: Convert Dovetail XXX.table.txt to AGP
## Usage: dovetail2agp.pl <-i XXX.table.txt>
# <-o organism (optional)>
# <-n assembly-name (optional)>
# <-d assembly-data (optional)>
# <-s assembly-source (optional)>
# <-D assembly-Description (optional)>
# <-h to display this help>";

####################
# declare variables
####################

my $inputfile = $opt_i || die $usage . "\n";
my $opt_organism = $opt_o || "NA";
my $opt_asmname = $opt_n || "NA";
my $opt_asmdate = $opt_d || "NA";
my $opt_asmsource = $opt_s || "NA";
my $opt_asmDescription = $opt_D || "NA";
defined($opt_h) && die $usage . "\n";

# open FILE, $inputfile or die $!;
my $outpath = dirname($inputfile);
my $outbase = basename( $inputfile, ".txt" );
my $FILE = OpenArchiveFile($inputfile) or die $!;

# declare variables
my $agpheader = join("\n",
	"##agp-version	2.0", 
	"# ORGANISM: ".$opt_organism, 
	"# ASSEMBLY NAME: ".$opt_asmname, 
	"# ASSEMBLY DATE: ".$opt_asmdate, 
	"# GENOME CENTER: ".$opt_asmsource, 
	"# DESCRIPTION: ".$opt_asmDescription
	);

# AGP format description
# col0: Obj_Name
# col1: Obj_Start
# col2: Obj_End 
# col3: PartNum 
# col4: Compnt_Type 'W' | 'N'
# col5: CompntId | Gap_Len
# col6: CompntStart | gap_type='scaffold'
# col7: CompntEnd | linkage='yes'
# col8: Orientation | 'map'

my $agpcolnames = join( "\t", 
	"# Obj_Name", 
	"Obj_Start", 
	"Obj_End", 
	"PartNum", 
	"Compnt_Type", 
	"CompntId|Gap_Len", 
	"CompntStart|gap_type", 
	"CompntEnd|linkage", 
	"Orientation|Linkage_evidence" );
	
my $count = 0;
my $first = 1;

################################
# load file into array of arrays
################################

my @data = ();
my %objCount = ();
my @obj = ();

while ( my $line = <$FILE> ) {
	# pass through header block
	if ( $line =~ /^#/ ) {
		next;
	}
	# split line in elements
	chomp($line);
	my @field = split("\t", $line);
	# reformat as AGP, add 1 to start coordinate
	my @agprow = ( $field[0], $field[5]+1, $field[6], "0", "W", $field[1], $field[2]+1, @field[3..4] );
	push @data, [ @agprow ];
	# store Object name and count
	$objCount{$field[0]}++;
	if ($objCount{$field[0]} == 1){
		push @obj, $field[0]
		};
	$count++;		
}	
print STDOUT "# ".$count." rows loaded from file, mapping to ".(keys %objCount)." Scaffolds\n";

# sort ObjName list naturally
my @sortedObj = natsort @obj;
our $cntgap = 0;

#################
# print out @agp
#################

# output AGP format
my $outfile = $outpath."/".$outbase.".agp";
open AGP, ">".$outfile || die $!;
print AGP $agpheader."\n";
print AGP $agpcolnames."\n";

# loop in @sortedObj and produce sorted blocks
foreach my $Obj_Name (@sortedObj) {
	# fetch all @data rows for that Obj_Name
	my @Objdata = grep {$_->[0] eq $Obj_Name} @data;
	# sort by coordinates
	my @sorted = arraysort(@Objdata);
	# add GAP rows and update partNum
	# extract first element from array
	my @first = shift @sorted;
	my @prev = @{$first[0]};
	my $PartNum = 1;
	# print first row to AGP
	$prev[3] = $PartNum;
	printagp(@prev);
	# loop through remaining rows
	foreach my $rec (@sorted) {
		my @field = @{$rec};
		# test if contiguous
		if ($field[1] == $prev[2]) {
			$PartNum ++;
			$field[3] = $PartNum;
			printagp(@field);
			@prev = @field;
		} else {
			# add GAP row
			$PartNum ++;
			my @gap = ( $field[0], $prev[2]+1, $field[1]-1, $PartNum );
			printgap(@gap);
			# also print $agp with PartNum incremented
			$PartNum++;			
			$field[3] = $PartNum;
			printagp(@field);
			@prev = @field;
		}
		#exit 0;
	}
}

# take care of handles neetly
undef $FILE;
close AGP;

print STDOUT "# AGP data of ".$cntgap." lines saved to ".$outfile."\n";

exit 0;

##############
#### Subs ####

sub printgap {
	my @a = @_;
	my @gap = ( @a[0..3], "N", $a[2]-$a[1]+1, "scaffold", "yes", "map" );
	print AGP join("\t", @gap)."\n";
	#print Dumper(@gap);
	$cntgap++;		
}

sub printagp {
	my @a = @_;
	print AGP join("\t", @a)."\n";
	#print Dumper(@a);		
	$cntgap++;
	}

sub arraysort {
	# sort the @results array by chromosome, start, end coordinate
	my @data = @_;
	my @sorted = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @data;
	return @sorted;
}

sub OpenArchiveFile {
	# $Filename passed in, handle to file passed out
	my $File = shift;	# filename
	my $FH;			  # file handle

	if ( $File =~ /.txt$/ ) {
		open( $FH, "cat $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /.txt.zip$/ ) {
		open( $FH, "unzip -p $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /(.txt.gzip|.bnx.gz)$/ ) {
		open( $FH, "gzip -dc $File | " ) or die("$!: can't open file $File");
	} else {
		die("$!: the file $File does seem to be a 'txt' file");
	}
	return $FH;
}
