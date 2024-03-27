#!/usr/bin/perl
#
# takes the genotype file and snp file and generates gemma input
# example: pattern_g_tcr_refugio_gs.txt snps_gs.txt
# example: color_g_tcr_refugio_gs.txt snps_gs.txt
#

$gf = shift(@ARGV);
$sf = shift(@ARGV);

open(IN, $sf) or die;
while(<IN>){
	chomp;
	push(@snps,$_);
}
close(IN);
open(IN, $gf);
$out = $gf;
$out =~ s/txt/geno/ or die "failed sub\n";
open(OUT, "> $out") or die;
while(<IN>){
	chomp;
	$s = shift(@snps);
	$s =~ s/ /_/;
	$o = "$s A T $_";
	$o =~ s/ /, /g;
	print OUT "$o\n";
}
close(IN);
close(OUT);
