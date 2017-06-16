#!/usr/bin/env perl

# BNGagpToLift - generate a lift file from a BioNano AGP file
# based on agp2Lift, see next line
# $Header: /projects/compbio/cvsroot/kent/src/utils/agpToLift,v 1.2 2007/03/09 19:29:00 angie Exp $
## first version: 2017-05-09
## extract original coordinate from AGP optical superscaffolds
## Bionano: deduce true coordinates from '_subseq_start:end'
#
# Stephane Plaisance (VIB-NC+BITS) 2017/05/09; v1.0
# visit our Git: https://github.com/Nucleomics-VIB

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

## INPUT formats (simplified)
# BioNano AGP format (GAP='N' in 5th column is excluded, only 'W' is used):
# object, object_beg, object_end, part_number, component_type (W),
#     component_id, component_beg, component_end, orientation (all should become '+')
# OUTPUT format compatible with liftUp
# offset oldName oldSize newName newSize

############################
# handle command parameters
############################

my $usage = "Aim: Convert BioNano AGP data to AGP format compatible with LiftUp (Kent tools). You must provide a AGP file with -i
## Usage: BNGagpToLift.pl <-i BioNano_hyrid-scaffoldAGP-file>
# <-h to display this help>";

getopts('i:h');
our ( $opt_i, $opt_h );
my $infile = $opt_i || die $usage."\n";
defined($opt_h) && die $usage."\n";

# open stream from AGP file
# open FILE, $inputfile or die $!;
my $FILE = OpenArchiveFile($infile) or die $!;
my $outpath = dirname($infile);

# remove possible suffixes from filename
my @sufx = ( ".agp", ".agp.gzip", ".agp.gz", ".agp.zip" );
my $outbase = basename( $infile, @sufx );

# create output handle
open LIFT, ">".$outpath."/".$outbase.".lift" || die $!;

# global ariables
my %chromSizes = ();
my @lines = ();
my $i;

# check for options

# read in AGP file
while ( my $line = <$FILE> ) {
	# ignore header lines
	next if $line =~ /^#/;

	# split line in elements
	chomp($line);
	my ( $chrom, $chromStart, $chromEnd, $ix, $type,
            $frag, $fragStart, $fragEnd, $strand ) = split("\t", $line);

	# skip N-gaps
	if ($type eq "N"){
		# do not process GAP lines
		next;
		}

    # data is sorted, overwrite 'chromSizes' to longest=last 'chromEnd'
    $chromSizes{$chrom} = $chromEnd;

	# debug
    # print STDERR $line,$chromEnd,"\n";

	# store in array for convertion
	$lines[$i++] = $line;
}

# close input stream
undef $FILE;

# process AGP lines, writing LIFT file, using chrom length info
for ($i=0; $i < @lines; $i++) {
    $_ = $lines[$i];
    my ( $chrom, $chromStart, $chromEnd, $ix, $type,
            $frag, $fragStart, $fragEnd, $strand ) = split /\s+/;

    # correct fragment name and coordinates
	if ($frag =~ m/(.*)(_subseq_)(.*)/g) {
		$frag = $1;
		my $coordinates = $3;
		my ( $subseqstart, $subseqend ) = ( split(/:/, $coordinates) );
		# update coordinates based on subseq info and strand
		$fragStart = $subseqstart-1;
		$fragEnd = $subseqend;
		}

    # output
    printf LIFT "%s\t%s\t%s\t%s\t%s\t%s\n",
            $chromStart - 1, $frag, ( $fragEnd - $fragStart + 1),
            $chrom, $chromSizes{$chrom}, $strand;
}

# take care of handles neetly
close LIFT;

exit 0;

##############
#### Subs ####

sub OpenArchiveFile {
	# $Filename passed in, handle to file passed out
	my $File = shift;	# filename
	my $FH;			  # file handle

	if ( $File =~ /.agp$/ ) {
		open( $FH, "cat $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /.agp.zip$/ ) {
		open( $FH, "unzip -p $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /(.agp.gzip|.agp.gz)$/ ) {
		open( $FH, "gzip -dc $File | " ) or die("$!: can't open file $File");
	} elsif ( $File =~ /(.agp.tgz|.agp.tar.gz)$/ ) {
		open( $FH, "tar -zxvf $File | " ) or die("$!: can't open file $File");
	} else {
		die("$!: the file $File does seem to be a 'agp' file");
	}
	return $FH;
}
