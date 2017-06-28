#!/usr/bin/perl

use Getopt::Std;
getopts "x:i:e:";

# https://github.com/tanghaibao/jcvi/issues/37
# https://github.com/tangerzhang/my_script/blob/master/bionano2Allmaps.pl

use Getopt::Std;
getopts "x:i:e:";


if ((!defined $opt_x) || (!defined $opt_i) || (!defined $opt_e)) {
    die "************************************************************************
    Usage: perl $0 -x xmap -i ref.fasta -e BspQI
      -h : help and usage.
      -x : xmap from hybrid assembly
      -i : ref.fasta
      -e : in silico digest enzyme, could be any dual enzyme
           e.g. -e BspQI,BbvCI
************************************************************************\n";
}else{
  print "************************************************************************\n";
  print "Version 1.0\n";
  print "Copyright to Tanger, tanger.zhang\@gmail.com\n";
  print "RUNNING...\n";
  print "************************************************************************\n";
	
	}
$ref_file = $opt_i;
$enzyme   = $opt_e;
$xmap     = $opt_x;
$enzyme   =~ s/,/ /g;

if($ref_file ne "ref.fasta"){
	$cmd = "ln -s $ref_file ./ref.fasta" ;
  system($cmd);
	}

# edit to your own path for fa2cmap_multi_color.pl
#$fa2cmap_multi_color='/Users/splaisan/git_repos/splaisan/BionanoThinkmate@NC/scripts/HybridScaffold/scripts/fa2cmap_multi_color.pl';
#$cmd = "perl $fa2cmap_multi_color -i ref.fasta -e $enzyme";

$fa2cmap='/Users/splaisan/git_repos/splaisan/BionanoThinkmate@NC/scripts/HybridScaffold/scripts/fa2cmap.pl';
$cmd = "perl $fa2cmap -i ref.fasta -n $enzyme";

system($cmd);

$enzyme =~ s/\s+/_/g;
$key_file = "fa2cmap/ref_".$enzyme."_0Kb_0labels_key.txt";

open(IN, $key_file) or die"";
while(<IN>){
	chomp;
	next if(/#/);
	next if(/CompntId/);
	@data = split(/\s+/,$_);
	$id   = $data[0];
	$scaf = $data[1];
	$iddb{$id} = $scaf;
	}
close IN;

open(IN, $xmap) or die"";
while(<IN>){
	chomp;
	next if(/#/);
	@data      = split(/\s+/,$_);
	$chr_id    = $data[1];
	$chr_posia = $data[3];
	$chr_posib = $data[4]; 
	$lg        = $data[2];
	$lg_posia  = $data[5];
	$lg_posib  = $data[6];
	$infordb{$lg}->{$lg_posia} = $chr_id."_".$chr_posia; 
	$infordb{$lg}->{$lg_posib} = $chr_id."_".$chr_posib; 
	}
close IN;
open(OUT, "> bionano.map.csv") or die"";
print OUT "Scaffold ID,scaffold position,LG,genetic position\n";
foreach $lg(sort {$a<=>$b} keys %infordb){
	foreach $lg_posi (sort {$a<=>$b} keys %{$infordb{$lg}}){
		($chr_id,$chr_posi) = split(/_/,$infordb{$lg}->{$lg_posi});
		$chrn = $iddb{$chr_id};
		print OUT "$chrn,$chr_posi,$lg,$lg_posi\n";
		}
	}
close OUT;